extern crate pretty_env_logger;

use clap::{Parser, ValueEnum};

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

    /// randomization strategy
    #[arg(value_enum, long, default_value_t = Randomizer::Shuffle)]
    pub random: Randomizer,

    /// overlap counting strategy
    #[arg(value_enum, long, default_value_t = Counter::All)]
    pub count: Counter,

    /// bed file of genome regions to mask (chrom<tab>start<tab>end)
    #[arg(long)]
    pub mask: Option<std::path::PathBuf>,

    /// randomize regions within each chromosome
    #[arg(long = "per-chrom", default_value_t = false)]
    pub per_chrom: bool,

    /// don't merge inputs' overlaps before processing
    #[arg(long = "no-merge-ovl", default_value_t = false)]
    pub no_merge: bool,

    /// do not swap A and B
    #[arg(long = "no-swap", default_value_t = false)]
    pub no_swap: bool,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Randomizer {
    // shuffle intervals allowing overlaps
    Shuffle = 0,
    // rotate intervals preserving order/spacing
    Circle = 1,
    // shuffle intervals without allowing overlaps
    Novl = 2,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Counter {
    // count number of overlaps
    All = 0,
    // count if any overlap
    Any = 1,
}

impl ArgParser {
    pub fn validate(&self) -> bool {
        let mut is_ok = true;
        if !self.bed_a.is_file() {
            error!("-A file doesn't exist");
            is_ok = false;
        }
        if !self.bed_b.is_file() {
            error!("-B file doesn't exist");
            is_ok = false;
        }
        if !self.genome.is_file() {
            error!("--genome file doesn't exist");
            is_ok = false;
        }

        if let Some(m) = &self.mask {
            if !m.is_file() {
                error!("--mask file doesn't exist");
                is_ok = false;
            }
        }

        if self.threads < 1 {
            warn!("need at least 1 thread");
            is_ok = false;
        }

        if (self.random == Randomizer::Novl) & self.no_merge {
            warn!("using `novl` without merged overlaps may cause errors");
            is_ok = false;
        }

        if self.num_times < 100 {
            warn!(
                "minimum p-value with {} is {}.",
                self.num_times,
                1.0 / ((self.num_times as f32) + 1.0)
            );
        }

        is_ok
    }
}
