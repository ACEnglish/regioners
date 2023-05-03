use clap::ValueEnum;
use serde::Serialize;
use rust_lapper::Lapper;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Overlapper {
    // count number of overlaps
    All,
    // count if any overlap
    Any,
}

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
