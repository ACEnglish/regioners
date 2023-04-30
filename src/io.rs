extern crate pretty_env_logger;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use rust_lapper::{Interval, Lapper};
pub type Iv = Interval<u64, u64>;

use std::collections::HashMap;
pub type MaskShift = HashMap<String, Lapper<u64, u64>>;

#[derive(Clone)]
pub struct GenomeShift {
    // start/end coordinates of a chromosome (for per-chrom) inside concat genome coordinates
    pub chrom: Lapper<u64, u64>,
    // chrom : how much to shift to get to concatenated genome coordinates
    pub shift: HashMap<String, u64>,
    // total span of the genome
    pub span: u64,
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_mask(file: &std::path::PathBuf) -> MaskShift {
    info!("parsing {}", file.display());
    let mut load: HashMap<String, Vec<Iv>> = HashMap::new();

    let mut num_mask = 0;

    if let Ok(lines) = read_lines(file) {
        for line in lines.flatten() {
            let collection: Vec<&str> = line.split('\t').collect();
            if collection.len() < 3 {
                error!("malformed bed line: {}", line);
                std::process::exit(1);
            }
            let chrom = collection[0].to_string();
            let m_start = collection[1].parse::<u64>().unwrap();
            let m_stop = collection[2].parse::<u64>().unwrap();

            if !load.contains_key(&chrom) {
                load.insert(chrom.clone(), Vec::<Iv>::new());
            }

            // load.get_mut(&chrom).map(|val| val.push(Iv{start: m_start, stop: m_stop, val: 0}));
            if let Some(val) = load.get_mut(&chrom) {
                val.push(Iv {
                    start: m_start,
                    stop: m_stop,
                    val: 0,
                })
            };
            num_mask += 1;
        }
    }

    let mut ret = MaskShift::new();
    for (key, val) in load.iter_mut() {
        val.sort();
        // can't assume they're sorted, so I gotta do this
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
    Read a two column file and turn into a GenomeShifter
    */
    info!("parsing {}", file.display());

    let mut load: Vec<Iv> = vec![];
    let mut m_shift: HashMap<String, u64> = HashMap::new();
    let mut cur_start: u64 = 0;
    let mut tot_masked: u64 = 0;

    if let Ok(lines) = read_lines(file) {
        for line in lines.flatten() {
            let collection: Vec<&str> = line.split('\t').collect();
            if collection.len() < 2 {
                error!("malformed genome line: {}", line);
                std::process::exit(1);
            }
            let chrom = collection[0].to_string();
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
            let size = collection[1].parse::<u64>().unwrap() - masked_bases;

            m_shift.insert(chrom.clone(), cur_start);
            load.push(Iv {
                start: cur_start,
                stop: cur_start + size,
                val: size,
            });
            cur_start += size;
            tot_masked += masked_bases;
        }
    }
    info!("loaded {} chromosomes", load.len());
    info!("total genome size: {}", cur_start);
    if tot_masked != 0 {
        info!("masked {} bases", tot_masked);
    }
    load.sort();

    GenomeShift {
        chrom: Lapper::new(load),
        shift: m_shift,
        span: cur_start,
    }
}

pub fn read_bed(
    file: &std::path::PathBuf,
    genome: &GenomeShift,
    mask: &Option<MaskShift>,
) -> Lapper<u64, u64> {
    /*
    Read a three column bed file and turn into a Lapper
    */
    info!("parsing {}", file.display());
    let mut ret: Vec<Iv> = vec![];
    let mut tot_size: u64 = 0;
    let mut num_masked = 0;
    let mut warned_chroms: Vec<String> = vec![];

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(file) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines.flatten() {
            let collection: Vec<&str> = line.split('\t').collect();
            if collection.len() < 3 {
                error!("malformed bed line: {}", line);
                std::process::exit(1);
            }
            let chrom = collection[0].to_string();

            if genome.shift.contains_key(&chrom) {
                let m_start = collection[1].parse::<u64>().unwrap();
                let m_stop = collection[2].parse::<u64>().unwrap();

                // if it intersects, get out
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

                // mask will r_shift
                ret.push(Iv {
                    start: m_start + r_shift - l_shift,
                    stop: m_stop + r_shift - l_shift,
                    val: 0,
                });
                tot_size += m_stop - m_start;
            } else if !warned_chroms.contains(&chrom) {
                // but only once
                warn!("{} missing from --genome and won't be loaded", chrom);
                warned_chroms.push(chrom);
            }
        }
    }
    info!("loaded {} intervals", ret.len());
    info!("masked {} intervals", num_masked);
    info!("total span: {}", tot_size);
    ret.sort();

    Lapper::new(ret)
}
