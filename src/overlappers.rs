//! Overlap counters
use clap::ValueEnum;
use rust_lapper::Lapper;
use serde::Serialize;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Overlapper {
    /// Count number of overlaps
    All,
    /// Count if any overlap
    Any,
}

/// For each interval in A, count any or all overlaps with B
impl Overlapper {
    pub fn ovl(&self, a_intv: &Lapper<u64, u64>, b_intv: &Lapper<u64, u64>) -> u64 {
        match self {
            /* Return number of A intervals intersecting a B intervals */
            Overlapper::All => a_intv
                .iter()
                .map(|i| match b_intv.find(i.start, i.stop).next() {
                    Some(_) => 1,
                    None => 0,
                })
                .sum(),
            /* Return number of b intervals intersecting each of a's intervals */
            Overlapper::Any => a_intv
                .iter()
                .map(|i| b_intv.find(i.start, i.stop).count() as u64)
                .sum(),
        }
    }
}
