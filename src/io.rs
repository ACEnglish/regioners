extern crate pretty_env_logger;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use rust_lapper::{Interval, Lapper};
pub type Iv = Interval<u64, u64>;

use std::collections::HashMap;
type GenomeShift = HashMap::<String, u64>;

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_bed(file: &std::path::PathBuf, genome: &GenomeShift) -> Lapper<u64, u64> {
    /*
    Read a three column bed file and turn into a Lapper
    */
    info!("parsing {}", file.display());
    let mut ret: Vec<Iv> = vec![];
    let mut tot_size:u64 = 0;

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(file) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(cur_line) = line {
                let collection: Vec<&str> = cur_line.split("\t").collect();
                if collection.len() < 3 {
                    error!("malformed bed line: {}", cur_line);
                    std::process::exit(1);
                }
                let chrom = collection[0].to_string();
                // raise error on key doesn't exist
                let m_shift = genome[&chrom];
                let m_start = collection[1].parse::<u64>().unwrap() + m_shift;
                let m_stop = collection[2].parse::<u64>().unwrap() + m_shift;
                ret.push(Iv{start: m_start,
                            stop: m_stop,
                            val: 0});
                tot_size += m_stop - m_start;
            }
        }
    }
    info!("loaded {} intervals", ret.len());
    info!("total span: {}", tot_size);
    return Lapper::new(ret);
}

pub fn read_genome(file: &std::path::PathBuf) -> (GenomeShift, u64) {
    /*
    Read a two column file and turn into a GenomeShifter
    */
    info!("parsing {}", file.display());

    let mut ret = GenomeShift::new();
    let mut cur_start:u64 = 0;

    if let Ok(lines) = read_lines(file) {
        for line in lines {
            if let Ok(cur_line) = line {
                let collection: Vec<&str> = cur_line.split("\t").collect();
                if collection.len() < 2 {
                    error!("malformed genome line: {}", cur_line);
                    std::process::exit(1);
                }
                let chrom = collection[0].to_string();
                let size = collection[1].parse::<u64>().unwrap();
                ret.insert(chrom, cur_start);
                cur_start += size;
            }
        }
    }
    info!("loaded {} chromosomes", ret.keys().len());
    info!("total genome size: {}", cur_start);
    return (ret, cur_start);
}
