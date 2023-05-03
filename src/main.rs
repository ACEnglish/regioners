extern crate pretty_env_logger;
#[macro_use]
extern crate log;

use std::fs::File;
use std::io::prelude::*;
use std::thread::{Builder, JoinHandle};

use clap::Parser;
use rust_lapper::Lapper;
use serde_json::json;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

mod cli;
mod gapbreaks;
mod io;
mod randomizers;

use crate::cli::ArgParser;
use crate::io::{read_bed, read_genome, read_mask};
use crate::randomizers::{circle_intervals, novl_intervals, shuffle_intervals};

// **********
// Helpers
fn all_overlap_count(a_intv: &Lapper<u64, u64>, b_intv: &Lapper<u64, u64>) -> u64 {
    /* Return number of b intervals intersecting each of a's intervals */
    a_intv.iter()
        .map(|i| b_intv.find(i.start, i.stop).count() as u64)
        .sum()
}

fn any_overlap_count(a_intv: &Lapper<u64, u64>, b_intv: &Lapper<u64, u64>) -> u64 {
    /* Return number of A intervals intersecting a B intervals */
    a_intv.iter()
        .map(|i| match b_intv.find(i.start, i.stop).next() {
            Some(_) => 1,
            None => 0,
        })
        .sum()
}

fn mean_std(v: &[u64]) -> (f64, f64) {
    /* Calculate mean and standard deviation */
    let n = v.len() as f64;
    let mean = v.iter().sum::<u64>() as f64 / n;
    let variance = v.iter().map(|x| (*x as f64 - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();

    (mean, std_dev)
}

fn count_permutations(o_count: u64, obs: &[u64], alt: char) -> f64 {
    /* Return number of permutations the observed count is (g)reater or (l)ess than */
    match alt {
        'l' => obs.iter().map(|i| (o_count >= *i) as u8 as f64).sum(),
        'g' => obs.iter().map(|i| (o_count <= *i) as u8 as f64).sum(),
        _ => panic!("Invalid alt"),
    }
}
// **********

fn main() -> std::io::Result<()> {
    pretty_env_logger::formatted_timed_builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    // IO
    let args = ArgParser::parse();
    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let mask = args.mask.map(|p| read_mask(&p));
    let mut genome = read_genome(&args.genome, &mask);
    let mut a_intv = read_bed(&args.bed_a, &genome, &mask);
    let mut b_intv = read_bed(&args.bed_b, &genome, &mask);

    // Setup
    if !args.no_merge {
        info!("merging overlaps");
        a_intv.merge_overlaps();
        b_intv.merge_overlaps();
    }
    let swapped = match !args.no_swap & (a_intv.len() > b_intv.len()) {
        true => {
            info!("swapping A for shorter B");
            std::mem::swap(&mut a_intv, &mut b_intv);
            true
        }
        false => false,
    };
    let overlapper = match args.count {
        cli::Counter::Any => any_overlap_count,
        cli::Counter::All => all_overlap_count,
    };
    let randomizer = match args.random {
        cli::Randomizer::Circle => circle_intervals,
        cli::Randomizer::Shuffle => shuffle_intervals,
        cli::Randomizer::Novl => {
            genome.make_gap_budget(&a_intv, &args.per_chrom);
            novl_intervals
        }
    };

    // profiling
    /*let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();*/

    // Processing
    let initial_overlap_count: u64 = overlapper(&a_intv, &b_intv);
    info!("observed : {}", initial_overlap_count);
    
    let progs = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .progress_chars("##-");
    
    let chunk_size: u32 = ((args.num_times as f32) / (args.threads as f32)).ceil() as u32;
    let mut pb : Vec<ProgressBar> = vec![];
    for _i in 0..args.threads {
        let p = progs.add(ProgressBar::new(chunk_size.into()));
        p.set_style(sty.clone());
        pb.push(p);
    }

    let handles: Vec<JoinHandle<Vec<u64>>> = (0..args.threads)
        .map(|i| {
            let m_a = a_intv.clone();
            let m_b = b_intv.clone();
            let m_g = genome.clone();
            let m_i = i.clone() as usize;
            let m_p = pb[m_i].clone();
            let builder = Builder::new().name(format!("regionRs-{}", i));
            // Send chunk to thread
            let start_iter = (i as u32) * chunk_size;
            let stop_iter = std::cmp::min(start_iter + chunk_size, args.num_times);
            builder
                .spawn(move || {
                    (start_iter..stop_iter).map(|_| { m_p.inc(1); overlapper(&randomizer(&m_a, &m_g, args.per_chrom), &m_b) }).collect()
                })
                .unwrap()
        })
        .collect();

    // Collect
    let mut all_counts: Vec<u64> = vec![];
    for handle in handles {
        let result: Vec<u64> = handle.join().unwrap();
        all_counts.extend(result);
    }
    progs.clear().unwrap();
    /*if let Ok(report) = guard.report().build() { println!("report: {:?}", &report); };*/

    // Calculate
    let (mu, sd) = mean_std(&all_counts);
    let alt = match (initial_overlap_count as f64) < mu {
        true => 'l',
        false => 'g',
    };
    let g_count = count_permutations(initial_overlap_count, &all_counts, alt);
    let p_val = (g_count + 1.0) / ((all_counts.len() as f64) + 1.0);
    let z_score = match (initial_overlap_count == 0) & (mu == 0.0) {
        true => {
            warn!("z_score cannot be computed");
            0.0
        }
        false => ((initial_overlap_count as f64) - mu) / sd,
    };

    // Output
    info!("perm mu: {}", mu);
    info!("perm sd: {}", sd);
    info!("alt hypo : {}", alt);
    info!("p-val : {}", p_val);

    let data = json!({"pval": p_val,
                      "zscore": z_score,
                      "obs":initial_overlap_count,
                      "perm_mu": mu,
                      "perm_sd": sd,
                      "alt": alt,
                      "n": args.num_times,
                      "swapped": swapped,
                      "no_merge": args.no_merge,
                      "random": args.random as u8,
                      "counter": args.count as u8,
                      "A_cnt" : a_intv.len(),
                      "B_cnt" : b_intv.len(),
                      "per_chrom": args.per_chrom,
                      "perms": all_counts});
    let json_str = serde_json::to_string(&data).unwrap();
    let mut file = File::create(args.output)?;
    file.write_all(json_str.as_bytes())
}
