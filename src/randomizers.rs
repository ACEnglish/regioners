//! Interval randomization techniques
use clap::ValueEnum;
use rust_lapper::Lapper;
use serde::Serialize;
use tinyrand::{Rand, RandRange, Seeded, StdRand};
use tinyrand_std::clock_seed::ClockSeed;

use crate::gapbreaks::GapBreaks;
use crate::io::{GenomeShift, Iv};

/// Holds the interval randomizations
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Randomizer {
    /// shuffle intervals allowing overlaps
    Shuffle,
    /// rotate intervals preserving order/spacing
    Circle,
    /// shuffle intervals without allowing overlaps
    Novl,
}

impl Randomizer {
    pub fn ize(
        &self,
        intv: &Lapper<u64, u64>,
        genome: &GenomeShift,
        per_chrom: bool,
    ) -> Lapper<u64, u64> {
        Lapper::<u64, u64>::new((match self {
            Randomizer::Circle => circle_intervals,
            Randomizer::Shuffle => shuffle_intervals,
            Randomizer::Novl => match genome.gap_budget {
                Some(_) => novl_intervals,
                None => panic!("Cannot run novl randomizer without gap_budget in genome"),
            },
        })(intv, genome, per_chrom))
    }
}

/// Randomly move each interval to a new position
fn shuffle_intervals(intv: &Lapper<u64, u64>, genome: &GenomeShift, per_chrom: bool) -> Vec<Iv> {
    let mut rand = StdRand::seed(ClockSeed::default().next_u64());
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
            let new_position = rand.next_range(lower..(upper - (i.stop - i.start)));
            Iv {
                start: new_position,
                stop: new_position + (i.stop - i.start),
                val: 0,
            }
        })
        .collect()
}

/// Randomly shift all intervals downstream with wrap-around
fn circle_intervals(intv: &Lapper<u64, u64>, genome: &GenomeShift, per_chrom: bool) -> Vec<Iv> {
    let mut rand = StdRand::seed(ClockSeed::default().next_u64());
    let mut ret: Vec<Iv> = Vec::<Iv>::with_capacity(intv.len());

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
    ret
}

/// Randomly move each interval to new position without overlapping them
fn novl_intervals(intv: &Lapper<u64, u64>, genome: &GenomeShift, per_chrom: bool) -> Vec<Iv> {
    let mut ret: Vec<Iv> = Vec::<Iv>::with_capacity(intv.len());

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
    ret
}

/// Shift each interval a set amount. Intervals shifted to < 0 are trimmed/removed
pub fn shift_intervals(intv: &Lapper<u64, u64>, shift: i64) -> Lapper<u64, u64> {
    Lapper::<u64, u64>::new(
        intv.iter()
            .map(|i| {
                let new_position = i.start as i64 + shift;
                let (new_start, new_end) = if new_position < 0 {
                    if (i.stop as i64) + shift <= 0 {
                        (0, 0)
                    } else {
                        (0, (i.stop - i.start) as i64 + new_position)
                    }
                } else {
                    (new_position, new_position + (i.stop - i.start) as i64)
                };
                Iv {
                    start: new_start as u64,
                    stop: new_end as u64,
                    val: 0,
                }
            })
            .collect(),
    )
}
