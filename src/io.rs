//! Input file parsers
extern crate pretty_env_logger;

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use rust_lapper::{Interval, Lapper};

pub type Iv = Interval<u64, u64>;
pub type MaskShift = HashMap<String, Lapper<u64, u64>>;

#[derive(Clone)]
pub struct GenomeShift {
    // start/end coordinates of a chromosome (for per-chrom) inside concat genome coordinates
    pub chrom: Lapper<u64, u64>,
    // chrom : how much to shift to get to concatenated genome coordinates
    pub shift: HashMap<String, u64>,
    // total span of the genome
    pub span: u64,
    pub gap_budget: Option<HashMap<u64, u64>>,
}

impl GenomeShift {
    pub fn make_gap_budget(&mut self, intervals: &Lapper<u64, u64>, per_chrom: &bool) {
        let mut ret = HashMap::<u64, u64>::new();
        match per_chrom {
            false => {
                ret.insert(0, self.span - intervals.cov());
            }
            true => {
                for i in self.chrom.iter() {
                    ret.insert(
                        i.start,
                        intervals
                            .find(i.start, i.stop)
                            .map(|p| p.stop - p.start)
                            .sum(),
                    );
                }
            }
        }
        self.gap_budget = Some(ret);
    }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
type FileHandler = io::Result<io::Lines<io::BufReader<File>>>;
pub fn read_lines<P>(filename: P) -> FileHandler
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

struct BedParser {
    /* Read tab delimited bed files while ensuring entries have start < end.
    It also ensures entries are sorted */
    file: std::path::PathBuf,
    prev_chrom: String,
    prev_start: u64,
}

impl BedParser {
    pub fn new(path: &Path) -> Self {
        Self {
            file: path.to_path_buf(),
            prev_chrom: String::new(),
            prev_start: 0,
        }
    }

    pub fn parse(&mut self, three_cols: bool) -> Vec<(String, u64, u64)> {
        if let Ok(lines) = read_lines(&mut self.file) {
            lines
                .flatten()
                .map(|line| {
                    let collection: Vec<&str> = line.split('\t').collect();
                    let n_cols = if three_cols { 3 } else { 2 };
                    if collection.len() < n_cols {
                        error!("malformed bed line: {}", line);
                        std::process::exit(1);
                    }
                    let chrom = collection[0].to_string();
                    let m_start = collection[1].parse::<u64>().unwrap();
                    let m_stop = if three_cols {
                        collection[2].parse::<u64>().unwrap()
                    } else {
                        m_start + 1
                    };

                    if chrom != self.prev_chrom {
                        self.prev_chrom = chrom.clone();
                        self.prev_start = 0;
                    }

                    if m_stop <= m_start {
                        error!("malformed bed line: stop <= start {}", line);
                        std::process::exit(1);
                    }
                    if m_start < self.prev_start {
                        error!(
                            "bed file unordered `sort -k3n -k1,2n` offending line {}",
                            line
                        );
                        std::process::exit(1);
                    }
                    (chrom, m_start, m_stop)
                })
                .collect()
        } else {
            panic!("unable to read bed file");
        }
    }
}

pub fn read_mask(file: &Path) -> MaskShift {
    /* read bed file into Mask Shift */
    info!("parsing {}", file.display());
    let mut load: HashMap<String, Vec<Iv>> = HashMap::new();
    let mut num_mask = 0;

    let mut m_parser = BedParser::new(file);
    for (chrom, m_start, m_stop) in m_parser.parse(true).into_iter() {
        if !load.contains_key(&chrom) {
            load.insert(chrom.clone(), Vec::<Iv>::new());
        }

        if let Some(val) = load.get_mut(&chrom) {
            val.push(Iv {
                start: m_start,
                stop: m_stop,
                val: 0,
            })
        };
        num_mask += 1;
    }

    let mut ret = MaskShift::new();
    for (key, val) in load.iter_mut() {
        ret.insert(key.clone(), Lapper::new(val.clone()));
    }

    info!(
        "loaded {} mask regions on {} chromosomes",
        num_mask,
        ret.keys().len()
    );

    ret
}

pub fn read_genome(file: &std::path::PathBuf, mask: &Option<MaskShift>) -> GenomeShift {
    /*
    Read a two column genome into a GenomeShifter
    */
    info!("parsing {}", file.display());

    let mut load: Vec<Iv> = vec![];
    let mut m_shift: HashMap<String, u64> = HashMap::new();
    let mut cur_start: u64 = 0;
    let mut tot_masked: u64 = 0;

    let mut m_parser = BedParser::new(file);
    for (chrom, mut size, _) in m_parser.parse(false) {
        let masked_bases = match mask {
            Some(m) => {
                if m.contains_key(&chrom) {
                    m[&chrom].cov()
                } else {
                    0
                }
            }
            None => 0,
        };
        size -= masked_bases;

        m_shift.insert(chrom.clone(), cur_start);
        load.push(Iv {
            start: cur_start,
            stop: cur_start + size,
            val: size,
        });
        cur_start += size;
        tot_masked += masked_bases;
    }
    info!("loaded {} chromosomes", load.len());
    info!("total genome size: {}", cur_start);
    if tot_masked != 0 {
        info!("masked {} bases", tot_masked);
    }

    GenomeShift {
        chrom: Lapper::new(load),
        shift: m_shift,
        span: cur_start,
        gap_budget: None,
    }
}

pub fn read_bed(file: &Path, genome: &GenomeShift, mask: &Option<MaskShift>) -> Lapper<u64, u64> {
    /*
    Read bed file into a Lapper
    */
    info!("parsing {}", file.display());
    let mut ret: Vec<Iv> = vec![];
    let mut tot_size: u64 = 0;
    let mut num_masked = 0;
    let mut warned_chroms: Vec<String> = vec![];

    let mut m_parser = BedParser::new(file);
    for (chrom, m_start, m_stop) in m_parser.parse(true).into_iter() {
        if !genome.shift.contains_key(&chrom) {
            // only warn once
            if !warned_chroms.contains(&chrom) {
                warn!("{} missing from --genome and won't be loaded", chrom);
                warned_chroms.push(chrom);
            }
            continue
        }
        let skip = match mask {
            Some(m) => {
                if m.contains_key(&chrom) {
                    m[&chrom].find(m_start, m_stop).count() != 0
                } else {
                    false
                }
            }
            None => false,
        };
        if skip {
            num_masked += 1;
            continue;
        }

        let r_shift = genome.shift[&chrom];
        let mut l_shift = 0;
        if mask.is_some() {
            let t_mask = mask.as_ref().expect("MASK SOME");
            if t_mask.contains_key(&chrom) {
                for i in t_mask[&chrom].find(0, m_start) {
                    l_shift += i.stop - i.start;
                }
            }
        }

        // end-to-end chrom coordinates will r_shift (increase).
        // maksed bases before start will l_shift (decrease).
        ret.push(Iv {
            start: m_start + r_shift - l_shift,
            stop: m_stop + r_shift - l_shift,
            val: 0,
        });
        tot_size += m_stop - m_start;
    }
    info!("loaded {} intervals", ret.len());
    info!("masked {} intervals", num_masked);
    info!("total span: {}", tot_size);

    Lapper::new(ret)
}
