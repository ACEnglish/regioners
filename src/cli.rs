extern crate pretty_env_logger;

use clap::Parser;

#[derive(Parser)]
// A rust implementation of regioneR for interval overlap permutation testing
#[command(author = "ACEnglish", version)]
pub struct ArgParser {
    /// chromosome lengths (chrom<tab>length)
    #[arg(short, long)]
    pub genome: std::path::PathBuf,

    /// bed file of regions (chrom<tab>start<tab>end)
    #[arg(short = 'A')]
    pub bed_a: std::path::PathBuf,

    /// bed file of regions (chrom<tab>start<tab>end)
    #[arg(short = 'B')]
    pub bed_b: std::path::PathBuf,

    /// number of permutations to perform
    #[arg(short, long = "num-times", default_value_t = 100)]
    pub num_times: u32,

    /// output json file
    #[arg(short, long)]
    pub output: std::path::PathBuf,

    /// number of threads to use
    #[arg(short, long, default_value_t = 1)]
    pub threads: u8,

    /// bed file of genome regions to mask (chrom<tab>start<tab>end)
    #[arg(long)]
    pub mask: Option<std::path::PathBuf>,

    /// do not swap A and B
    #[arg(long = "no-swap", default_value_t = false)]
    pub no_swap: bool,

    /// count any overlap instead of number of overlaps
    #[arg(long, default_value_t = false)]
    pub any: bool,

    /* Depricated for Day 2
    /// don't allow overlapping entries during randomization
    #[arg(long = "no-overlap", default_value_t = false)]
    pub no_overlaps: bool,
    */
    /// randomize regions within each chromosome
    #[arg(long = "per-chrom", default_value_t = false)]
    pub per_chrom: bool,

    /// use circularization for randomization
    #[arg(long, default_value_t = false)]
    pub circle: bool,
    /* Depricated for Day 2
    /// maximum randomization retries before fail
    #[arg(long = "max-retry", default_value_t = 250)]
    pub max_retry: u64
    */
}

pub fn validate_args(args: &ArgParser) -> bool {
    let mut is_ok = true;
    if !args.bed_a.is_file() {
        error!("-A file doesn't exist");
        is_ok = false;
    }
    if !args.bed_b.is_file() {
        error!("-B file doesn't exist");
        is_ok = false;
    }
    if !args.genome.is_file() {
        error!("--genome file doesn't exist");
        is_ok = false;
    }
    /*if !args.mask.is_file() { an optional so, if it is there I have to check it..
        error!("--mask file doesn't exist");
        is_ok = false;
    }*/
    if args.num_times < 100 {
        warn!(
            "minimum p-value with {} is {}.",
            args.num_times,
            1.0 / ((args.num_times as f32) + 1.0)
        );
    }
    // TODO: thread checks, mask check

    is_ok
}
