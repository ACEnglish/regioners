use plotters::prelude::*;
use serde_json::Value;
use std::error::Error;
use std::f64::consts::PI;
use std::path::Path;

const W: u32 = 1800;
const H: u32 = 900;
const FONT: &str = "sans-serif";
const FONT_TICK: u32 = 22;
const FONT_LABEL: u32 = 26;
const FONT_TITLE: u32 = 30;

// ── KDE helpers ────────────────────────────────────────────────────────────

fn scotts_bandwidth(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    let mean = data.iter().sum::<f64>() / n;
    let std = (data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt();
    1.06 * std * n.powf(-0.2)
}

fn kde(data: &[f64], bw: f64, x: f64) -> f64 {
    let n = data.len() as f64;
    data.iter()
        .map(|xi| (-(x - xi).powi(2) / (2.0 * bw * bw)).exp())
        .sum::<f64>()
        / (n * bw * (2.0 * PI).sqrt())
}

/// Returns (bin_left, bin_right, density) for each bin, normalised to density.
fn density_bins(data: &[f64], n_bins: usize) -> Vec<(f64, f64, f64)> {
    let lo = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let width = (hi - lo) / n_bins as f64;
    let mut counts = vec![0usize; n_bins];
    for &x in data {
        let i = ((x - lo) / width).floor() as usize;
        counts[i.min(n_bins - 1)] += 1;
    }
    let total = data.len() as f64;
    (0..n_bins)
        .map(|i| {
            let l = lo + i as f64 * width;
            (l, l + width, counts[i] as f64 / (total * width))
        })
        .collect()
}

// ── Public entry point ──────────────────────────────────────────────────────

/// Derives SVG paths from the JSON output path and writes both plots.
pub fn plot_results(data: &Value, json_path: &Path) -> Result<(), Box<dyn Error>> {
    let stem = json_path.file_stem().unwrap().to_string_lossy();
    let dir  = json_path.parent().unwrap_or(Path::new("."));
    plot_perm_dist(data, &dir.join(format!("{stem}_perms.png")))?;
    plot_local_z(data,   &dir.join(format!("{stem}_localz.png")))?;
    Ok(())
}

// ── Plot 1: permutation distribution ───────────────────────────────────────

fn plot_perm_dist(data: &Value, path: &Path) -> Result<(), Box<dyn Error>> {
    let test = &data["test"];
    let perms: Vec<f64> = test["perms"]
        .as_array().ok_or("missing 'perms'")?
        .iter().filter_map(|v| v.as_f64()).collect();
    let obs = test["observed"].as_f64().ok_or("missing 'observed'")?;

    // x range — always include the observed value
    let x_lo = perms.iter().cloned().fold(f64::INFINITY,     f64::min).min(obs);
    let x_hi = perms.iter().cloned().fold(f64::NEG_INFINITY, f64::max).max(obs);
    let pad  = (x_hi - x_lo) * 0.05;
    let (x_lo, x_hi) = (x_lo - pad, x_hi + pad);

    let bins = density_bins(&perms, 30);
    let bw   = scotts_bandwidth(&perms);
    let kde_pts: Vec<(f64, f64)> = (0..=400)
        .map(|i| {
            let x = x_lo + (x_hi - x_lo) * i as f64 / 400.0;
            (x, kde(&perms, bw, x))
        })
        .collect();

    // y ceiling — max of histogram and KDE, with headroom for the annotation
    let y_max = bins.iter().map(|b| b.2)
        .chain(kde_pts.iter().map(|p| p.1))
        .fold(0.0_f64, f64::max)
        * 1.15;

    let root = BitMapBackend::new(path, (W, H)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(30)
        .x_label_area_size(70)
        .y_label_area_size(90)
        .build_cartesian_2d(x_lo..x_hi, 0.0..y_max)?;

    chart.configure_mesh()
        .x_desc("Intersection Count")
        .y_desc("Permutation Density")
        .axis_desc_style((FONT, FONT_LABEL).into_font())
        .label_style((FONT, FONT_TICK).into_font())
        .draw()?;

    // Gray histogram bars
    chart.draw_series(bins.iter().map(|&(l, r, d)| {
        Rectangle::new(
            [(l, 0.0), (r, d)],
            ShapeStyle { color: RGBAColor(170, 170, 170, 0.8), filled: true, stroke_width: 1 },
        )
    }))?;

    // Black KDE curve
    chart.draw_series(LineSeries::new(kde_pts, &BLACK))?;

    // Blue vertical line at observed
    chart.draw_series(LineSeries::new(
        [(obs, 0.0), (obs, y_max * 0.95)],
        ShapeStyle::from(&BLUE).stroke_width(2),
    ))?;

    // Annotation: wheat box + rotated label, centred on y_max/2
    let y_mid   = y_max / 2.0;
    let box_w   = (x_hi - x_lo) * 0.05;
    let box_h   = y_max * 0.30;
    /* chart.draw_series(std::iter::once(Rectangle::new(
        [(obs, y_mid - box_h / 2.0), (obs + box_w, y_mid + box_h / 2.0)],
        ShapeStyle { color: RGBAColor(245, 222, 179, 0.9), filled: true, stroke_width: 1 },
    )))?; */

    let label = "observed intersections";

    // Plotters has no text-measure API, so estimate pixel width from char count.
    // Sans-serif characters are roughly 0.55× as wide as they are tall.
    // Then convert pixels → data units using the known axis/canvas dimensions.
    let plot_h_px = (H - 2 * 70 - 90) as f64; // H - 2*margin - x_label_area_size
    let text_px   = label.len() as f64 * FONT_TICK as f64 * 0.55;
    let text_span = text_px * (y_max / plot_h_px);

    chart.draw_series(std::iter::once(Text::new(
        label,
        (obs + box_w * 0.15, y_mid + box_h - text_span / 2.0),
        (FONT, FONT_TICK)
            .into_font()
            .transform(FontTransform::Rotate90)
            .color(&BLACK),
    )))?;

    root.present()?;
    info!("Intersection plot saved");
    Ok(())
}

// ── Plot 2: local z-scores ─────────────────────────────────────────────────

fn plot_local_z(data: &Value, path: &Path) -> Result<(), Box<dyn Error>> {
    let lz     = &data["localZ"];
    let window = lz["window"].as_i64().ok_or("missing 'window'")? as i32;
    let step   = lz["step"].as_i64().ok_or("missing 'step'")? as i32;
    let shifts: Vec<f64> = lz["shifts"]
        .as_array().ok_or("missing 'shifts'")?
        .iter().filter_map(|v| v.as_f64()).collect();

    let xs: Vec<i32> = (-window..window).step_by(step as usize).collect();
    let n = xs.len().min(shifts.len());
    let points: Vec<(i32, f64)> =
        xs[..n].iter().copied().zip(shifts[..n].iter().copied()).collect();

    let x_lo  = *xs.first().unwrap();
    let x_hi  = *xs.last().unwrap();
    let y_lo  = shifts[..n].iter().cloned().fold(f64::INFINITY,     f64::min);
    let y_hi  = shifts[..n].iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_pad = (y_hi - y_lo) * 0.1;

    let root = BitMapBackend::new(path, (W, H)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Local z-score values", ("sans-serif", 22))
        .margin(30)
        .x_label_area_size(70)
        .y_label_area_size(90)
        .build_cartesian_2d(x_lo..x_hi, (y_lo - y_pad)..(y_hi + y_pad))?;

    chart.configure_mesh()
        .x_desc("Shift")
        .y_desc("z-score")
        .axis_desc_style((FONT, FONT_LABEL).into_font())
        .label_style((FONT, FONT_TICK).into_font())
        .draw()?;

    chart.draw_series(LineSeries::new(points, &BLUE))?;

    root.present()?;
    info!("Local Z plot saved");
    Ok(())
}
