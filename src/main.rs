use anyhow::{bail, Context, Result};
use clap::Parser;
use gpx::read;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::PathBuf;

#[derive(clap::ValueEnum, Clone, Debug)]
enum PinShape {
    Circle,
    Teardrop,
}

#[derive(Parser, Debug)]
#[command(name = "gpx2svg", about = "Convert GPX tracks/routes to an SVG outline")]
struct Args {
    /// Input GPX file
    #[arg(short, long)]
    input: PathBuf,

    /// Output SVG file (defaults to stdout if omitted)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// SVG canvas width in pixels
    #[arg(long, default_value_t = 1200)]
    width: u32,

    /// SVG canvas height in pixels
    #[arg(long, default_value_t = 1200)]
    height: u32,

    /// Padding in pixels around the path
    #[arg(long, default_value_t = 40.0)]
    padding: f64,

    /// Stroke color (any valid CSS/SVG color, e.g. "#ff00aa" or "black")
    #[arg(long, default_value = "#000000")]
    stroke: String,

    /// Stroke width in pixels
    #[arg(long, default_value_t = 2.0)]
    stroke_width: f64,

    /// Fill color for the path (usually "none" for outlines)
    #[arg(long, default_value = "none")]
    fill: String,

    /// If set, close the path (useful for loop tracks)
    #[arg(long, default_value_t = false)]
    close_path: bool,

    /// Simplify geometry with RDP using epsilon in *pixels* (0 disables)
    #[arg(long, default_value_t = 0.0)]
    simplify_px: f64,

    /// Draw a start marker at the first point
    #[arg(long, default_value_t = false)]
    start_pin: bool,

    /// Draw an end marker at the last point
    #[arg(long, default_value_t = false)]
    end_pin: bool,

    /// Marker radius in pixels
    #[arg(long, default_value_t = 6.0)]
    pin_radius: f64,

    /// Start marker fill color
    #[arg(long, default_value = "#16a34a")] // green-ish
    start_pin_fill: String,

    /// End marker fill color
    #[arg(long, default_value = "#dc2626")] // red-ish
    end_pin_fill: String,

    /// Marker stroke color
    #[arg(long, default_value = "#ffffff")]
    pin_stroke: String,

    /// Marker stroke width
    #[arg(long, default_value_t = 2.0)]
    pin_stroke_width: f64,
    /// Shape of the start/end pins
    #[arg(long, value_enum, default_value = "circle")]
    pin_shape: PinShape,
}

#[derive(Clone, Copy, Debug)]
struct Pt {
    x: f64,
    y: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let file = File::open(&args.input)
        .with_context(|| format!("Failed to open input GPX: {}", args.input.display()))?;
    let reader = BufReader::new(file);
    let gpx = read(reader).context("Failed to parse GPX")?;

    // Collect points from tracks and routes.
    // geo_types::Point uses x=lon, y=lat.
    let mut latlons: Vec<(f64, f64)> = Vec::new(); // (lat, lon)

    for trk in &gpx.tracks {
        for seg in &trk.segments {
            for pt in &seg.points {
                let p = pt.point();
                latlons.push((p.y(), p.x())); // (lat, lon)
            }
        }
    }

    for rte in &gpx.routes {
        for pt in &rte.points {
            let p = pt.point();
            latlons.push((p.y(), p.x())); // (lat, lon)
        }
    }

    if latlons.len() < 2 {
        bail!("No track/route points found (need at least 2 points).");
    }

    // --- Project lat/lon to a flat plane (equirectangular) ---
    let mean_lat = latlons.iter().map(|(lat, _)| lat).sum::<f64>() / latlons.len() as f64;
    let cos_lat = mean_lat.to_radians().cos();

    let pts: Vec<Pt> = latlons
        .into_iter()
        .map(|(lat, lon)| Pt {
            x: lon * cos_lat,
            y: lat,
        })
        .collect();

    // --- Fit transform (projected -> pixel space) ---
    let (min_x, max_x, min_y, max_y) = bounds(&pts);

    let w = args.width as f64;
    let h = args.height as f64;

    let inner_w = (w - 2.0 * args.padding).max(1.0);
    let inner_h = (h - 2.0 * args.padding).max(1.0);

    let span_x = (max_x - min_x).abs().max(1e-12);
    let span_y = (max_y - min_y).abs().max(1e-12);

    let scale = (inner_w / span_x).min(inner_h / span_y);

    let geom_w = span_x * scale;
    let geom_h = span_y * scale;
    let offset_x = args.padding + (inner_w - geom_w) / 2.0;
    let offset_y = args.padding + (inner_h - geom_h) / 2.0;

    // Map projected points to pixel-space SVG coords (invert Y).
    let mut px_pts: Vec<Pt> = pts
        .iter()
        .map(|p| Pt {
            x: offset_x + (p.x - min_x) * scale,
            y: offset_y + (max_y - p.y) * scale, // invert Y
        })
        .collect();

    // Capture marker positions BEFORE simplification (RDP keeps endpoints anyway,
    // but this avoids any future changes affecting marker selection).
    let start_marker = px_pts.first().copied();
    let end_marker = px_pts.last().copied();

    // --- Simplify in pixel space ---
    if args.simplify_px > 0.0 {
        px_pts = rdp_simplify(&px_pts, args.simplify_px);
    }

    let d = to_svg_path_d(&px_pts, args.close_path);

    // Optional marker SVG
    let mut markers = String::new();
if args.start_pin == true {
    if let Some(p) = start_marker {
        let marker = match args.pin_shape {
            PinShape::Circle => {
                svg_circle(
                            p,
                            args.pin_radius,
                            &args.start_pin_fill,
                            &args.pin_stroke,
                            args.pin_stroke_width,
                        )
            },
            PinShape::Teardrop => {
                svg_teardrop(
                                p,
                                args.pin_radius,
                                &args.start_pin_fill,
                                &args.pin_stroke,
                                args.pin_stroke_width,
                            )
            },
        };
        markers.push_str(&marker);
        markers.push('\n');
    }
}
    if args.end_pin {
        if let Some(p) = end_marker {
            markers.push_str(&svg_circle(
                p,
                args.pin_radius,
                &args.end_pin_fill,
                &args.pin_stroke,
                args.pin_stroke_width,
            ));
            markers.push('\n');
        }
    }

    let svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w0}" height="{h0}" viewBox="0 0 {w0} {h0}">
  <path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" stroke-linejoin="round" stroke-linecap="round" />
  {markers}</svg>
"#,
        w0 = args.width,
        h0 = args.height,
        d = d,
        fill = escape_attr(&args.fill),
        stroke = escape_attr(&args.stroke),
        sw = args.stroke_width,
        markers = markers
    );

    // Save to file if --output is provided; otherwise print to stdout.
    match args.output {
        Some(path) => {
            let mut f = File::create(&path)
                .with_context(|| format!("Failed to create output: {}", path.display()))?;
            f.write_all(svg.as_bytes())?;
        }
        None => {
            print!("{svg}");
        }
    }

    Ok(())
}

fn bounds(pts: &[Pt]) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;

    for p in pts {
        min_x = min_x.min(p.x);
        max_x = max_x.max(p.x);
        min_y = min_y.min(p.y);
        max_y = max_y.max(p.y);
    }
    (min_x, max_x, min_y, max_y)
}

fn to_svg_path_d(pts: &[Pt], close_path: bool) -> String {
    let mut s = String::new();
    if let Some(first) = pts.first() {
        s.push_str(&format!("M {:.2} {:.2}", first.x, first.y));
        for p in &pts[1..] {
            s.push_str(&format!(" L {:.2} {:.2}", p.x, p.y));
        }
        if close_path {
            s.push_str(" Z");
        }
    }
    s
}

fn svg_circle(p: Pt, r: f64, fill: &str, stroke: &str, stroke_width: f64) -> String {
    format!(
        r#"<circle cx="{:.2}" cy="{:.2}" r="{:.2}" fill="{}" stroke="{}" stroke-width="{}" />"#,
        p.x,
        p.y,
        r,
        escape_attr(fill),
        escape_attr(stroke),
        stroke_width
    )
}
fn svg_teardrop(
    p: Pt,
    size: f64,
    fill: &str,
    stroke: &str,
    stroke_width: f64,
) -> String {
    // Teardrop pointing down, anchored so the "tip" is at (p.x, p.y)
    let r = size;
    let cy = p.y - r; // center of the circular part

    format!(
        r#"<path d="
M {x:.2} {tip_y:.2}
C {x1:.2} {y1:.2}, {x2:.2} {y1:.2}, {x2:.2} {cy:.2}
A {r:.2} {r:.2} 0 1 1 {x0:.2} {cy:.2}
C {x0:.2} {y1:.2}, {x1:.2} {y1:.2}, {x:.2} {tip_y:.2}
Z"
fill="{fill}" stroke="{stroke}" stroke-width="{sw}" />
"#,
        x = p.x,
        tip_y = p.y,
        cy = cy,
        r = r,
        x0 = p.x - r,
        x1 = p.x,
        x2 = p.x + r,
        y1 = cy + r * 1.4,
        fill = escape_attr(fill),
        stroke = escape_attr(stroke),
        sw = stroke_width
    )
}


fn escape_attr(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

// ---------------------
// RDP simplification
// ---------------------

fn rdp_simplify(points: &[Pt], epsilon: f64) -> Vec<Pt> {
    if points.len() <= 2 || epsilon <= 0.0 {
        return points.to_vec();
    }

    let mut keep = vec![false; points.len()];
    keep[0] = true;
    keep[points.len() - 1] = true;

    rdp_mark(points, 0, points.len() - 1, epsilon, &mut keep);

    points
        .iter()
        .zip(keep.iter())
        .filter_map(|(p, k)| if *k { Some(*p) } else { None })
        .collect()
}

fn rdp_mark(points: &[Pt], start: usize, end: usize, epsilon: f64, keep: &mut [bool]) {
    if end <= start + 1 {
        return;
    }

    let a = points[start];
    let b = points[end];

    let mut max_dist = -1.0_f64;
    let mut idx = 0usize;

    for i in (start + 1)..end {
        let d = perp_dist(points[i], a, b);
        if d > max_dist {
            max_dist = d;
            idx = i;
        }
    }

    if max_dist > epsilon {
        keep[idx] = true;
        rdp_mark(points, start, idx, epsilon, keep);
        rdp_mark(points, idx, end, epsilon, keep);
    }
}

fn perp_dist(p: Pt, a: Pt, b: Pt) -> f64 {
    let vx = b.x - a.x;
    let vy = b.y - a.y;
    let wx = p.x - a.x;
    let wy = p.y - a.y;

    let vv = vx * vx + vy * vy;
    if vv <= 1e-24 {
        return ((p.x - a.x).powi(2) + (p.y - a.y).powi(2)).sqrt();
    }

    let t = (wx * vx + wy * vy) / vv;
    let t = t.clamp(0.0, 1.0);

    let proj_x = a.x + t * vx;
    let proj_y = a.y + t * vy;

    ((p.x - proj_x).powi(2) + (p.y - proj_y).powi(2)).sqrt()
}

