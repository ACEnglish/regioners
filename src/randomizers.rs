use tinyrand::{Rand, RandRange, Seeded, StdRand};
use tinyrand_std::clock_seed::ClockSeed;

use crate::gapbreaks::GapBreaks;
use crate::io::{GenomeShift, Iv};
use rust_lapper::Lapper;

// ***********
// Randomizers
// ***********
pub fn shuffle_intervals(
    intv: &Lapper<u64, u64>,
    genome: &GenomeShift,
    per_chrom: bool,
) -> Lapper<u64, u64> {
    /*
        Randomly move each interval to new position
    */
    let mut rand = StdRand::seed(ClockSeed::default().next_u64());
    Lapper::<u64, u64>::new(
        intv.iter()
            .map(|i| {
                let (lower, upper) = if per_chrom {
                    match genome.chrom.find(i.start, i.stop).next() {
                        Some(b) => (b.start, b.stop),
                        None => panic!("Interval @ ({}, {}) not hitting genome", i.start, i.stop),
                    }
                } else {
                    (0, genome.span)
                };
                let shift = rand.next_range(lower..(upper - (i.stop - i.start)));
                Iv {
                    start: i.start + shift,
                    stop: i.stop + shift,
                    val: 0,
                }
            })
            .collect(),
    )
}

pub fn circle_intervals(
    intv: &Lapper<u64, u64>,
    genome: &GenomeShift,
    per_chrom: bool,
) -> Lapper<u64, u64> {
    /*
        Randomly shift all intervals downstream with wrap-around
    */
    let mut rand = StdRand::seed(ClockSeed::default().next_u64());
    let mut ret = Vec::<Iv>::new();

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
            ret.push(Iv {
                start: new_start - upper,
                stop: new_end - upper,
                val: 0,
            });
        } else if new_end > upper {
            ret.extend(vec![
                Iv {
                    start: new_start,
                    stop: upper,
                    val: 0,
                },
                Iv {
                    start: lower,
                    stop: new_end - upper,
                    val: 0,
                },
            ]);
        } else {
            ret.push(Iv {
                start: new_start,
                stop: new_end,
                val: 0,
            });
        }
    }
    Lapper::<u64, u64>::new(ret)
}

pub fn novl_intervals(
    intv: &Lapper<u64, u64>,
    genome: &GenomeShift,
    per_chrom: bool,
) -> Lapper<u64, u64> {
    /*
        Randomly move each interval to new position without overlapping them
    */
    let mut ret: Vec<Iv> = vec![];

    let spans = match per_chrom {
        true => genome.chrom.intervals.clone(),
        false => vec![Iv {
            start: 0,
            stop: genome.span,
            val: 0,
        }],
    };

    for subi in spans {
        let m_gap: u64 = match &genome.gap_budget {
            Some(g) => g[&subi.start],
            None => panic!("How are you using the gap_budget without making it first?"),
        };

        let mut cur_intervals: Vec<(bool, u64)> = GapBreaks::new(m_gap).collect();
        cur_intervals.extend(
            intv.find(subi.start, subi.stop)
                .map(|i| (true, i.stop - i.start)),
        );
        fastrand::shuffle(&mut cur_intervals);

        let mut cur_pos = subi.start;
        for i in cur_intervals {
            if i.0 {
                ret.push(Iv {
                    start: cur_pos,
                    stop: cur_pos + i.1,
                    val: 0,
                });
            }
            cur_pos += i.1;
        }
    }
    Lapper::<u64, u64>::new(ret)
}
