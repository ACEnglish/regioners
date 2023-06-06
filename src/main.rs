//! A rust implementation of regioneR for interval overlap permutation testing
extern crate pretty_env_logger;
#[macro_use]
extern crate log;

use std::fs::File;
use std::io::prelude::*;
use std::sync::Arc;
use std::thread::JoinHandle;

use clap::Parser;
#[cfg(feature = "progbars")]
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use serde_json::json;

mod cli;
mod gapbreaks;
mod io;
mod overlappers;
mod randomizers;
mod stats;

use crate::cli::ArgParser;
use crate::io::{read_bed, read_genome, read_mask};
use crate::randomizers::Randomizer;
use crate::stats::{LocalZscore, PermTest};

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

    let mask = args.mask.as_ref().map(|p| read_mask(p));
    let mut genome = read_genome(&args.genome, &mask);
    let mut a_intv = read_bed(&args.bed_a, &genome, &mask);
    let mut b_intv = read_bed(&args.bed_b, &genome, &mask);

    // Setup
    if !args.no_merge {
        info!("merging overlaps");
        a_intv.merge_overlaps();
        b_intv.merge_overlaps();
    }
    let a_count = a_intv.len();
    let b_count = b_intv.len();
    let swapped = if !args.no_swap & (a_count > b_count) {
        info!("swapping A for shorter B");
        std::mem::swap(&mut a_intv, &mut b_intv);
        true
    } else {
        false
    };
    if args.random == Randomizer::Novl {
        genome.make_gap_budget(&a_intv, &args.per_chrom)
    }
    // Won't need to change again. Can pass pointers to threads
    let genome = Arc::new(genome);
    let a_intv = Arc::new(a_intv);
    let b_intv = Arc::new(b_intv);

    // profiling
    /*let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();*/

    // Processing
    let initial_overlap_count: u64 = args.count.ovl(&a_intv, &b_intv);
    info!("observed : {}", initial_overlap_count);

    #[cfg(feature = "progbars")]
    let (progs, pb) = {
        let chunk_size: u32 = ((args.num_times as f32) / (args.threads as f32)).ceil() as u32;
        let progs = MultiProgress::new();
        let sty = ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-");

        let mut pb: Vec<ProgressBar> = vec![];
        for _i in 0..args.threads {
            let p = progs.add(ProgressBar::new(chunk_size.into()));
            p.set_style(sty.clone());
            pb.push(p);
        }
        (progs, pb)
    };

    let handles: Vec<JoinHandle<Vec<u64>>> = (0..args.threads)
        .map(|i| {
            let m_a = a_intv.clone();
            let m_b = b_intv.clone();
            let m_g = genome.clone();
            #[cfg(feature = "progbars")]
            let m_p = pb[i as usize].clone();

            std::thread::spawn(move || {
                ((i as u32)..args.num_times)
                    .step_by(args.threads as usize)
                    .map(|_| {
                        #[cfg(feature = "progbars")]
                        m_p.inc(1);
                        args.count.ovl(&args.random.ize(&m_a, &m_g, args.per_chrom), &m_b)
                    })
                    .collect()
            })
        })
        .collect();

    // Collect
    let mut perm_counts = Vec::<u64>::with_capacity(args.num_times as usize);
    for handle in handles {
        perm_counts.extend(handle.join().unwrap());
    }
    #[cfg(feature = "progbars")]
    progs.clear().unwrap();
    /*if let Ok(report) = guard.report().build() { println!("report: {:?}", &report); };*/

    // Calculate
    let test = PermTest::new(initial_overlap_count, perm_counts);
    let local_zscores = LocalZscore::new(&a_intv, &b_intv, args.count, args.window, args.step, &test);

    // Output
    info!("perm mu: {}", test.mean);
    info!("perm sd: {}", test.std_dev);
    info!("alt hypo : {}", test.alt);
    info!("p-val : {}", test.p_val);
    let data = json!({"test": test,
                      "swapped": swapped,
                      "no_merge": args.no_merge,
                      "random": args.random,
                      "count": args.count,
                      "A_cnt" : a_count,
                      "B_cnt" : b_count,
                      "per_chrom": args.per_chrom,
                      "localZ": local_zscores,
    });
    let mut file = File::create(args.output)?;
    file.write_all(serde_json::to_string(&data).unwrap().as_bytes())
}
