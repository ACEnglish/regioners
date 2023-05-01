extern crate pretty_env_logger;
#[macro_use]
extern crate log;

use std::fs::File;
use std::io::prelude::*;
use std::thread;

use tinyrand::{Rand, RandRange, Seeded, StdRand};
use tinyrand_std::clock_seed::ClockSeed;

use clap::Parser;
use rust_lapper::Lapper;
use serde_json::json;

mod cli;
mod io;

// ***********
// Randomizers
// ***********
fn shuffle_intervals(
    intv: &Lapper<u64, u64>,
    genome: &io::GenomeShift,
    per_chrom: bool,
) -> Lapper<u64, u64> {
    /*
        Randomly move each interval to new position
    */
    let mut ret = Vec::<io::Iv>::new();
    let mut rand = StdRand::seed(ClockSeed::default().next_u64());

    for i in intv.iter() {
        let (lower, upper) = if per_chrom {
            match genome.chrom.find(i.start, i.stop).next() {
                Some(b) => (b.start, b.stop),
                None => panic!("Interval @ ({}, {}) not hitting genome", i.start, i.stop),
            }
        } else {
            (0, genome.span)
        };
        let span = i.stop - i.start;
        let shift = rand.next_range(lower..(upper - span));
        ret.push(io::Iv {
            start: i.start + shift,
            stop: i.stop + shift,
            val: 0,
        });
    }
    ret.sort();

    Lapper::<u64, u64>::new(ret)
}

fn circle_intervals(
    intv: &Lapper<u64, u64>,
    genome: &io::GenomeShift,
    per_chrom: bool,
) -> Lapper<u64, u64> {
    /*
        Randomly shift all intervals downstream with wrap-around
    */
    let mut rand = StdRand::seed(ClockSeed::default().next_u64());
    let mut ret = Vec::<io::Iv>::new();

    let genome_shift: u64 = rand.next_range(0..(genome.span));

    for i in intv.iter() {
        let (lower, upper, shift) = if per_chrom {
            match genome.chrom.find(i.start, i.stop).next() {
                Some(b) => (b.start, b.stop, genome_shift % b.val),
                None => panic!("Interval @ ({}, {}) not hitting genome", i.start, i.stop),
            }
        } else {
            (0, genome.span, genome_shift)
        };

        let new_start: u64 = i.start + shift;
        let new_end: u64 = i.stop + shift;

        if new_start >= upper {
            ret.push(io::Iv {
                start: new_start - upper,
                stop: new_end - upper,
                val: 0,
            });
        } else if new_end > upper {
            ret.push(io::Iv {
                start: new_start,
                stop: upper,
                val: 0,
            });
            ret.push(io::Iv {
                start: lower,
                stop: new_end - upper,
                val: 0,
            });
        } else {
            ret.push(io::Iv {
                start: new_start,
                stop: new_end,
                val: 0,
            });
        }
    }
    ret.sort();

    Lapper::<u64, u64>::new(ret)
}

// **********
// Overlapers
// **********
fn get_num_overlap_count(a_lap: &Lapper<u64, u64>, b_lap: &Lapper<u64, u64>) -> u64 {
    /*
        Return number of b intervals intersecting each of a's intervals
    */
    let mut inter_cnt: u64 = 0;
    for i in a_lap.iter() {
        inter_cnt += b_lap.find(i.start, i.stop).count() as u64;
    }

    inter_cnt
}

fn get_any_overlap_count(a_lap: &Lapper<u64, u64>, b_lap: &Lapper<u64, u64>) -> u64 {
    /*
        Return number of a intervals intersecting b intervals
    */
    let mut inter_cnt: u64 = 0;
    for i in a_lap.iter() {
        if b_lap.find(i.start, i.stop).next().is_some() {
            inter_cnt += 1;
        }
    }

    inter_cnt
}

// *******
// Helpers
// *******
fn mean_std(v: &[u64]) -> (f64, f64) {
    let n = v.len() as f64;
    let mean = v.iter().sum::<u64>() as f64 / n;

    let variance = v.iter().map(|x| (*x as f64 - mean).powi(2)).sum::<f64>() / n;

    let std_dev = variance.sqrt();

    (mean, std_dev)
}

fn count_permutations(o_count: u64, obs: &Vec<u64>, alt: char) -> f64 {
    /*
        Return number of permutations the observed count is (g)reater or (l)ess than
    */
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

    g_count
}

fn main() -> std::io::Result<()> {
    // ToDo - Debug?
    pretty_env_logger::formatted_timed_builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = cli::ArgParser::parse();
    if !cli::validate_args(&args) {
        error!("please fix arguments");
        std::process::exit(1);
    }

    // This was a manual implementation of map `match args.mask { Some(p) => Some(io::read_mask(&p)), None => None }`
    let mask = args.mask.map(|p| io::read_mask(&p));

    let genome = io::read_genome(&args.genome, &mask);

    let mut a_lapper = io::read_bed(&args.bed_a, &genome, &mask);
    let mut b_lapper = io::read_bed(&args.bed_b, &genome, &mask);
    if args.merge_overlaps {
        info!("merging overlaps");
        a_lapper.merge_overlaps();
        b_lapper.merge_overlaps();
    }

    let mut was_swapped = false;
    if !args.no_swap & (a_lapper.len() > b_lapper.len()) {
        info!("swapping A for shorter B");
        std::mem::swap(&mut a_lapper, &mut b_lapper);
        was_swapped = true;
    }

    let overlapper = if args.any {
        get_any_overlap_count
    } else {
        get_num_overlap_count
    };
    let initial_overlap_count: u64 = overlapper(&a_lapper, &b_lapper);
    info!("{} intersections", initial_overlap_count);

    let mut handles = Vec::new();
    let chunk_size: u32 = ((args.num_times as f32) / (args.threads as f32)).ceil() as u32;

    for i in 0..args.threads {
        let m_a = a_lapper.clone();
        let m_b = b_lapper.clone();
        let m_genome = genome.clone();
        //let m_ovl_arg = args.no_overlaps.clone();
        //let m_mxr_arg = args.max_retry.clone();

        let start_iter = (i as u32) * chunk_size;
        let stop_iter = std::cmp::min(start_iter + chunk_size, args.num_times);
        let handle = thread::spawn(move || {
            let mut m_counts: Vec<u64> = vec![];
            //let input_intv_cnt = m_a.len();
            for _j in start_iter..stop_iter {
                //let mut is_placed = false;
                //let mut num_attempts = 0;
                //while num_attempts < m_mxr_arg {
                    let new_intv = if args.circle {
                        circle_intervals(&m_a, &m_genome, args.per_chrom)
                    } else {
                        shuffle_intervals(&m_a, &m_genome, args.per_chrom)
                    };

                    //if args.circle | !m_ovl_arg | (new_intv.len() == input_intv_cnt) {
                        //is_placed = true;
                        m_counts.push(overlapper(&new_intv, &m_b));
                        //break;
                    //}
                    //num_attempts += 1;
                //}

                //if !is_placed {
                    //warn!("shuffling failed to produce --no-overlaps");
                //}
            }
            m_counts
        });
        handles.push(handle);
    }

    let mut all_counts: Vec<u64> = vec![];
    for handle in handles {
        let result = handle.join().unwrap();
        all_counts.extend(result);
    }

    let (mu, sd) = mean_std(&all_counts);
    info!("perm mu: {}", mu);
    info!("perm sd: {}", sd);

    let alt = if (initial_overlap_count as f64) < mu {
        'l'
    } else {
        'g'
    };
    info!("alt : {}", alt);
    let g_count = count_permutations(initial_overlap_count, &all_counts, alt);
    let p_val = (g_count + 1.0) / ((all_counts.len() as f64) + 1.0);
    info!("p-val : {}", p_val);

    let mut z_score = 0.0;
    if (initial_overlap_count == 0) & (mu == 0.0) {
        warn!("z_score cannot be computed");
    } else {
        z_score = ((initial_overlap_count as f64) - mu) / sd;
    }

    let data = json!({"pval": p_val,
                      "zscore": z_score,
                      "obs":initial_overlap_count,
                      "perm_mu": mu,
                      "perm_sd": sd,
                      "alt": alt,
                      "n": args.num_times,
                      "swapped": was_swapped,
                      "merged": args.merge_overlaps,
                      "A_cnt" : a_lapper.len(),
                      "B_cnt" : b_lapper.len(),
                      "any": args.any,
                      "per_chrom": args.per_chrom,
                      "circle": args.circle,
                      "perms": all_counts});
    let json_str = serde_json::to_string(&data).unwrap();

    // Write the JSON string to a file
    let mut file = File::create(args.output)?;
    file.write_all(json_str.as_bytes())
}
