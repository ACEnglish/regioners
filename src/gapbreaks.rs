//! Helper for breaking gaps during novl randomization
use tinyrand::{Rand, RandRange, Seeded, StdRand};
use tinyrand_std::ClockSeed;

/// When performing novl randomization, we break the uncovered spans of
/// the genome into pieces and shuffle them along with the intervals.
/// Truly random novl would break all the uncovered spans into 1bp pieces.
/// However, this is extremly inefficient. Instead, we break it into pieces
/// randomly between 1bp and some maximum size. If this maximum size was
/// the length of all uncovered spans, we could end up making one giant
/// uncovered span and then grouping all the intervals end-to-end.
/// By setting the maximum size to some percent of the uncovered span
/// we lower the bias towards making large gaps.
/// NOVLMAGIC is calculated as `int(1 / max_size_percent)`
const NOVLMAGIC: u64 = 10000;

/// Holds the length of uncovered spans in a genome and a private random
/// number generator. Can be iterated to generate random sized gaps.
/// Gaps are returned as `(false, size)` tuples so they be differentiated
/// from covered intervals (`(true, size)`) in [`randomizers::novl_intervals`]
///
/// Example:
///  
/// ```
/// let mut gapsizes: Vec<(bool, u64)> = GapBreaks::new(100).collect();
/// ```
pub struct GapBreaks {
    total_gap_size: u64,
    rand: StdRand,
}

impl GapBreaks {
    pub fn new(total_gap_size: u64) -> Self {
        Self {
            total_gap_size: total_gap_size,
            rand: StdRand::seed(ClockSeed::default().next_u64()),
        }
    }
}

impl Iterator for GapBreaks {
    type Item = (bool, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.total_gap_size > 0 {
            let g_l = self
                .rand
                .next_range(1..std::cmp::max(2, self.total_gap_size / NOVLMAGIC));
            self.total_gap_size -= g_l;
            return Some((false, g_l));
        }
        None
    }
}
