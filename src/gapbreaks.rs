/* Helpers for performing breaking gaps during novl randomization */
use tinyrand::{Rand, RandRange, Seeded, StdRand};
use tinyrand_std::ClockSeed;

const NOVLMAGIC: u64 = 10000;
/* When performing novl randomization, we break the uncovered spans of
*  the genome into pieces and shuffle them along with the intervals.
*  Truly random novl would break all the uncovered spans into 1bp pieces.
*  However, this is extremly inefficient. Instead, we break it into pieces
*  randomly between 1bp and some maximum size. If this maximum size was
*  the length of all uncovered spans, we could end up making one giant
*  uncovered span and then grouping all the intervals end-to-end.
*  By setting the maximum size to some percent of the uncovered span
*  we lower the bias towards making large gaps.
*  NOVLMAGIC is calculated as int(1 / max_size_percent)
*/

pub struct GapBreaks {
    budget: u64,
    rand: StdRand,
}

impl GapBreaks {
    pub fn new(m_budget: u64) -> Self {
        Self {
            budget: m_budget,
            rand: StdRand::seed(ClockSeed::default().next_u64()),
        }
    }
}

impl Iterator for GapBreaks {
    type Item = (bool, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.budget > 0 {
            let g_l = self
                .rand
                .next_range(1..std::cmp::max(2, self.budget / NOVLMAGIC));
            self.budget -= g_l;
            return Some((false, g_l));
        }
        None
    }
}
