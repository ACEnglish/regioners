use crate::randomizers::shift_intervals;
use crate::ArgParser;
use rust_lapper::Lapper;
use serde::Serialize;

/// Creates and holds permutation test results
#[derive(Serialize)]
pub struct PermTest {
    pub observed: u64,
    pub num_perms: f64,
    pub mean: f64,
    pub std_dev: f64,
    pub p_val: f64,
    pub z_score: f64,
    pub alt: char,
    pub perms: Vec<u64>,
}

impl PermTest {
    pub fn new(observed: u64, perms: Vec<u64>) -> Self {
        let n = perms.len() as f64;
        let mean = perms.iter().sum::<u64>() as f64 / n;
        let variance = perms
            .iter()
            .map(|x| (*x as f64 - mean).powi(2))
            .sum::<f64>()
            / n;
        let std_dev = variance.sqrt();
        let (alt, p_count): (char, f64) = if (observed as f64) < mean {
            (
                'l',
                perms.iter().map(|i| (*i < observed) as u8 as f64).sum(),
            )
        } else {
            (
                'g',
                perms.iter().map(|i| (*i > observed) as u8 as f64).sum(),
            )
        };
        let p_val = (p_count + 1.0) / (n + 1.0);
        let z_score = if (observed == 0) & (mean == 0.0) {
            warn!("z_score cannot be computed");
            0.0
        } else {
            ((observed as f64) - mean) / std_dev
        };
        PermTest {
            observed,
            p_val,
            num_perms: n,
            z_score,
            mean,
            std_dev,
            alt,
            perms,
        }
    }
}

/// Creates and holds local z-score results
#[derive(Serialize)]
pub struct LocalZscore {
    pub shifts: Vec<f64>,
    pub window: i64,
    pub step: u64,
}

impl LocalZscore {
    pub fn new(
        a_intv: &Lapper<u64, u64>,
        b_intv: &Lapper<u64, u64>,
        args: &ArgParser,
        test: &PermTest,
    ) -> Self {
        let shifts: Vec<f64> = (-args.window..args.window)
            .step_by(args.step as usize)
            .map(|i| {
                let observed = args.count.ovl(&shift_intervals(a_intv, i), b_intv);
                if (observed == 0) & (test.mean == 0.0) {
                    0.0
                } else {
                    ((observed as f64) - test.mean) / test.std_dev
                }
            })
            .collect();
        LocalZscore {
            shifts,
            window: args.window,
            step: args.step,
        }
    }
}
