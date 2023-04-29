extern crate pretty_env_logger;
#[macro_use] extern crate log;

use std::thread;
use std::fs::File;
use std::io::prelude::*;

use rand::Rng;
use clap::Parser;
use rust_lapper::{Lapper};
use serde_json::{json};

mod cli;
mod io;


fn mean_std(v: &[u64]) -> (f64, f64) {
    let n = v.len() as f64;
    let mean = v.iter().sum::<u64>() as f64 / n;

    let variance = v.iter()
        .map(|x| (*x as f64 - mean).powi(2))
        .sum::<f64>() / n;

    let std_dev = variance.sqrt();

    (mean, std_dev)
}

fn get_overlap_count(a_lap: &Lapper<u64, u64>, b_lap: &Lapper<u64, u64>) -> u64 {
    let mut inter_cnt: u64 = 0;
    let mut cursor = 0;
    for i in a_lap.iter() {
        inter_cnt += b_lap.seek(i.start, i.stop, &mut cursor).count() as u64;
    }
    return inter_cnt;
}

fn shuffle_intervals(intv: &Lapper<u64, u64>, genome_size: u64) -> Lapper<u64, u64> {
    let mut ret: Vec<io::Iv> = vec![];
    let mut rng = rand::thread_rng();

    for i in intv.iter() {
        let mut is_placed = false;
        let span = i.stop - i.start;
        let max_end = genome_size - span;
        while !is_placed { // Should maybe have a max-tries so we don't hang here forever...
            let new_pos: u64 = rng.gen_range(0..max_end);
            let new_end: u64 = new_pos + span;
            if intv.find(new_pos, new_pos + span).count() == 0 {
                ret.push(io::Iv{start: new_pos, stop: new_end, val: 0});
                is_placed = true;
            }
        }
    }
    ret.sort();
    return Lapper::new(ret);
}

fn count_permutations(o_count: u64, obs: &Vec<u64>, alt: char) -> f64 {
    let mut g_count: f64 = 0.0;
    if alt == 'l' {
        for i in obs {
            g_count += if o_count >= *i { 1.0 } else { 0.0 };
        }
    } else if alt == 'g' {
        for i in obs {
            g_count += if o_count <= *i { 1.0 } else { 0.0 };
        }
    }
    return g_count
}

fn main() {
    // ToDo - Debug?
    pretty_env_logger::formatted_timed_builder().filter_level(log::LevelFilter::Info).init();

    let args = cli::ArgParser::parse();
    if !cli::validate_args(&args) {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let (genome, genome_size) = io::read_genome(&args.genome);

    let a_lapper = io::read_bed(&args.bed_a, &genome);
    let b_lapper = io::read_bed(&args.bed_b, &genome);

    if a_lapper.len() > b_lapper.len() {
        info!("swapping A for shorter B");
        let (a_lapper, b_lapper) = (b_lapper.clone(), a_lapper.clone());
    }
    
    let initial_overlap_count = get_overlap_count(&a_lapper, &b_lapper);
    info!("{} intersections", initial_overlap_count);
    
    let mut handles = Vec::new();
    let chunk_size:u32 = ((args.num_times as f32) / (args.threads as f32)).ceil() as u32;
    
    for i in 0..args.threads {
        let m_alap_copy = a_lapper.clone();
        let m_blap_copy = b_lapper.clone();
        let m_gs_copy = genome_size.clone();
        // There's got to be a better way
        let start_iter = (i as u32) * chunk_size;
        let stop_iter = std::cmp::min(start_iter + chunk_size, args.num_times);
        let m_range = start_iter..stop_iter;
        let handle = thread::spawn(move || {
            let mut m_counts:Vec<u64> = vec![];
            for _j in m_range {
                let new_intv = shuffle_intervals(&m_alap_copy, m_gs_copy);
                m_counts.push(get_overlap_count(&new_intv, &m_blap_copy));
            }
            m_counts
        });
        handles.push(handle);
    }

    let mut all_counts:Vec<u64> = vec![];
    for handle in handles {
        let result = handle.join().unwrap();
        all_counts.extend(result);
    }

    let (mu, sd) = mean_std(&all_counts);
    info!("perm mu: {}", mu);
    info!("perm sd: {}", sd);

    let alt = if (initial_overlap_count as f64) < mu {'l'} else {'g'};
    let g_count = count_permutations(initial_overlap_count, &all_counts, alt);
    let p_val = (g_count + 1.0) / ((all_counts.len() as f64) + 1.0);
    info!("p-val : {}", p_val);

    let data = json!({"pval": p_val, 
                      "zscore": ((initial_overlap_count as f64) - mu) / sd,
                      "obs":initial_overlap_count,
                      "perm_mu": mu,
                      "perm_sd": sd,
                      "alt": alt,
                      "n": args.num_times,
                      "perms": all_counts});
    let json_str = serde_json::to_string(&data).unwrap();

    // Write the JSON string to a file
    let file = File::create(args.output);
    file.expect("REASON").write_all(json_str.as_bytes());
}
