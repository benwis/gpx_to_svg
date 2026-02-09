use anyhow::{bail, Context, Result};
use clap::{ArgAction, Parser, ValueEnum};
use gpx::read;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::PathBuf;

#[derive(ValueEnum, Clone, Debug)]
enum PinShape {
    Circle,
    Teardrop,
}

#[derive(Parser, Debug)]
#[command(name = "gpx2svg", about = "Convert GPX tracks/routes to an SVG outline + optional elevation profile")]
struct Args {
    /// Input GPX file
    #[arg(short, long)]
    input: PathBuf,

    /// Output route SVG file (defaults to stdout if omitted)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// SVG canvas width in pixels
    #[arg(long, default_value_t = 1200)]
    width: u32,

    /// SVG canvas height in pixels
    #[arg(long, default_value_t = 1200)]
    height: u32,

    /// Padding in pixels around the route drawing area
    #[arg(long, default_value_t = 40.0)]
    padding: f64,

    /// Route stroke color (CSS/SVG color)
    #[arg(long, default_value = "#000000")]
    stroke: String,

    /// Route stroke width in pixels
    #[arg(long, default_value_t = 2.0)]
    stroke_width: f64,

    /// Route fill (usually "none")
    #[arg(long, default_value = "none")]
    fill: String,

    /// Close the route path (useful for loops)
    #[arg(long, default_value_t = false)]
    close_path: bool,

    /// Simplify route geometry (epsilon in pixels; 0 disables)
    #[arg(long, default_value_t = 0.0)]
    simplify_px: f64,

    /// Draw a start marker at the first point
    #[arg(long, default_value_t = false)]
    start_pin: bool,

    /// Draw an end marker at the last point
    #[arg(long, default_value_t = false)]
    end_pin: bool,

    /// Shape of start/end pins
    #[arg(long, value_enum, default_value = "circle")]
    pin_shape: PinShape,

    /// Pin size (circle radius, teardrop size)
    #[arg(long, default_value_t = 6.0)]
    pin_radius: f64,

    /// Start pin fill color
    #[arg(long, default_value = "#16a34a")]
    start_pin_fill: String,

    /// End pin fill color
    #[arg(long, default_value = "#dc2626")]
    end_pin_fill: String,

    /// Pin stroke color
    #[arg(long, default_value = "#ffffff")]
    pin_stroke: String,

    /// Pin stroke width
    #[arg(long, default_value_t = 2.0)]
    pin_stroke_width: f64,

    // -------------------------
    // Elevation options
    // -------------------------

    /// Write elevation profile to a separate SVG file (optional)
    #[arg(long)]
    elevation_svg: Option<PathBuf>,

    /// Elevation SVG width (only used for --elevation-svg)
    #[arg(long, default_value_t = 1200)]
    elevation_width: u32,

    /// Elevation SVG height (only used for --elevation-svg)
    #[arg(long, default_value_t = 320)]
    elevation_height: u32,

    /// Elevation padding inside the elevation chart
    #[arg(long, default_value_t = 30.0)]
    elevation_padding: f64,

    /// Elevation line stroke color
    #[arg(long, default_value = "#111827")]
    elevation_stroke: String,

    /// Elevation line stroke width
    #[arg(long, default_value_t = 2.0)]
    elevation_stroke_width: f64,

    /// Elevation fill under the curve ("none" disables fill)
    #[arg(long, default_value = "none")]
    elevation_fill: String,

    /// Elevation axis/label color
    #[arg(long, default_value = "#374151")]
    elevation_label_color: String,

    /// Draw elevation axes (true/false)
    #[arg(long, action = ArgAction::Set, default_value_t = true)]
    elevation_axes: bool,

    /// Draw elevation labels (true/false)
    #[arg(long, action = ArgAction::Set, default_value_t = true)]
    elevation_labels: bool,

    /// Embed elevation chart into the same route SVG (bottom section)
    #[arg(long, default_value_t = false)]
    embed_elevation: bool,

    /// Fraction of total height used by the route when embedding elevation (0.0–1.0)
    #[arg(long, default_value_t = 0.7)]
    embed_route_ratio: f64,

    /// When embedding elevation: match the elevation chart width to the route outline width
    #[arg(long, default_value_t = true, action = ArgAction::Set)]
    elevation_match_route_width: bool,
}

#[derive(Clone, Copy, Debug)]
struct Pt {
    x: f64,
    y: f64,
}

#[derive(Clone, Copy, Debug)]
struct Sample {
    lat: f64,
    lon: f64,
    ele: Option<f64>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // --- Read GPX ---
    let file = File::open(&args.input)
        .with_context(|| format!("Failed to open input GPX: {}", args.input.display()))?;
    let reader = BufReader::new(file);
    let gpx = read(reader).context("Failed to parse GPX")?;

    // --- Collect samples from tracks + routes ---
    // geo_types::Point: x=lon, y=lat
    let mut samples: Vec<Sample> = Vec::new();

    for trk in &gpx.tracks {
        for seg in &trk.segments {
            for pt in &seg.points {
                let p = pt.point();
                samples.push(Sample {
                    lat: p.y(),
                    lon: p.x(),
                    ele: pt.elevation,
                });
            }
        }
    }
    for rte in &gpx.routes {
        for pt in &rte.points {
            let p = pt.point();
            samples.push(Sample {
                lat: p.y(),
                lon: p.x(),
                ele: pt.elevation,
            });
        }
    }

    if samples.len() < 2 {
        bail!("No track/route points found (need at least 2 points).");
    }

    // -------------------------
    // Compute embedded layout
    // -------------------------
    let total_w = args.width as f64;
    let total_h = args.height as f64;

    let (route_view_h, elev_view_y, elev_view_h) = if args.embed_elevation {
        let r = args.embed_route_ratio.clamp(0.1, 0.9);
        let route_h = total_h * r;
        let elev_h = total_h - route_h;
        (route_h, route_h, elev_h)
    } else {
        (total_h, total_h, 0.0)
    };

    // -------------------------
    // ROUTE SVG (outline)
    // -------------------------

    // Equirectangular projection; use mean latitude to reduce distortion.
    let mean_lat = samples.iter().map(|s| s.lat).sum::<f64>() / samples.len() as f64;
    let cos_lat = mean_lat.to_radians().cos();

    let projected: Vec<Pt> = samples
        .iter()
        .map(|s| Pt {
            x: s.lon * cos_lat,
            y: s.lat,
        })
        .collect();

    // Fit route into the route viewport: [0..width] x [0..route_view_h]
    let (min_x, max_x, min_y, max_y) = bounds(&projected);

    let inner_w = (total_w - 2.0 * args.padding).max(1.0);
    let inner_h = (route_view_h - 2.0 * args.padding).max(1.0);

    let span_x = (max_x - min_x).abs().max(1e-12);
    let span_y = (max_y - min_y).abs().max(1e-12);

    // Uniform scale preserve aspect ratio
    let scale = (inner_w / span_x).min(inner_h / span_y);

    // Center route within route viewport
    let geom_w = span_x * scale;
    let geom_h = span_y * scale;
    let offset_x = args.padding + (inner_w - geom_w) / 2.0;
    let offset_y = args.padding + (inner_h - geom_h) / 2.0;

    // Project -> pixel space (invert Y)
    let mut px_pts: Vec<Pt> = projected
        .iter()
        .map(|p| Pt {
            x: offset_x + (p.x - min_x) * scale,
            y: offset_y + (max_y - p.y) * scale,
        })
        .collect();

    // Capture marker positions (in pixel space)
    let start_marker = px_pts.first().copied();
    let end_marker = px_pts.last().copied();

    // Compute route pixel span (for matching elevation width)
    let (route_left, route_right) = px_x_span(&px_pts);

    // Simplify in pixel space (if requested)
    if args.simplify_px > 0.0 {
        px_pts = rdp_simplify(&px_pts, args.simplify_px);
    }

    let path_d = to_svg_path_d(&px_pts, args.close_path);

    // Markers
    let mut markers = String::new();
    if args.start_pin {
        if let Some(p) = start_marker {
            markers.push_str(&match args.pin_shape {
                PinShape::Circle => svg_circle(
                    p,
                    args.pin_radius,
                    &args.start_pin_fill,
                    &args.pin_stroke,
                    args.pin_stroke_width,
                ),
                PinShape::Teardrop => svg_teardrop(
                    p,
                    args.pin_radius,
                    &args.start_pin_fill,
                    &args.pin_stroke,
                    args.pin_stroke_width,
                ),
            });
            markers.push('\n');
        }
    }
    if args.end_pin {
        if let Some(p) = end_marker {
            markers.push_str(&match args.pin_shape {
                PinShape::Circle => svg_circle(
                    p,
                    args.pin_radius,
                    &args.end_pin_fill,
                    &args.pin_stroke,
                    args.pin_stroke_width,
                ),
                PinShape::Teardrop => svg_teardrop(
                    p,
                    args.pin_radius,
                    &args.end_pin_fill,
                    &args.pin_stroke,
                    args.pin_stroke_width,
                ),
            });
            markers.push('\n');
        }
    }

    // Embedded elevation fragment (optional)
    let mut elevation_fragment = String::new();
    if args.embed_elevation {
        let span_override = if args.elevation_match_route_width {
            Some((route_left, route_right))
        } else {
            None
        };

        // If not enough elevation data, skip embedding gracefully
        if let Ok(frag) = build_elevation_svg_fragment(
            &samples,
            0.0,
            elev_view_y,
            total_w,
            elev_view_h,
            args.elevation_padding,
            &args.elevation_stroke,
            args.elevation_stroke_width,
            &args.elevation_fill,
            &args.elevation_label_color,
            args.elevation_axes,
            args.elevation_labels,
            span_override,
        ) {
            elevation_fragment = frag;
        }
    }

    // Final route SVG
    let route_svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">
  <g class="route">
    <path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" stroke-linejoin="round" stroke-linecap="round" />
    {markers}
  </g>
  {elev}
</svg>
"#,
        w = args.width,
        h = args.height,
        d = path_d,
        fill = escape_attr(&args.fill),
        stroke = escape_attr(&args.stroke),
        sw = args.stroke_width,
        markers = markers,
        elev = elevation_fragment
    );

    // Write route SVG
    match &args.output {
        Some(path) => {
            let mut f = File::create(path)
                .with_context(|| format!("Failed to create output: {}", path.display()))?;
            f.write_all(route_svg.as_bytes())?;
        }
        None => {
            print!("{route_svg}");
        }
    }

    // Optional separate elevation SVG
    if let Some(path) = &args.elevation_svg {
        let elev_svg = build_elevation_svg_full(
            &samples,
            args.elevation_width as f64,
            args.elevation_height as f64,
            args.elevation_padding,
            &args.elevation_stroke,
            args.elevation_stroke_width,
            &args.elevation_fill,
            &args.elevation_label_color,
            args.elevation_axes,
            args.elevation_labels,
        )?;

        let mut f = File::create(path)
            .with_context(|| format!("Failed to create elevation output: {}", path.display()))?;
        f.write_all(elev_svg.as_bytes())?;
    }

    Ok(())
}

// -------------------------
// Elevation rendering
// -------------------------

fn build_elevation_svg_full(
    samples: &[Sample],
    w: f64,
    h: f64,
    pad: f64,
    stroke: &str,
    sw: f64,
    fill: &str,
    label_color: &str,
    draw_axes: bool,
    draw_labels: bool,
) -> Result<String> {
    let frag = build_elevation_svg_fragment(
        samples,
        0.0,
        0.0,
        w,
        h,
        pad,
        stroke,
        sw,
        fill,
        label_color,
        draw_axes,
        draw_labels,
        None, // standalone chart: no route span to match
    )?;

    Ok(format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w0}" height="{h0}" viewBox="0 0 {w0} {h0}">
{frag}
</svg>
"#,
        w0 = w as u32,
        h0 = h as u32,
        frag = frag
    ))
}

fn build_elevation_svg_fragment(
    samples: &[Sample],
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    pad: f64,
    stroke: &str,
    sw: f64,
    fill: &str,
    label_color: &str,
    draw_axes: bool,
    draw_labels: bool,
    plot_left_right: Option<(f64, f64)>,
) -> Result<String> {
    // Build (distance_m, elevation_m) series, skipping points without elevation.
    let mut series: Vec<(f64, f64)> = Vec::new();
    let mut dist_m = 0.0;

    if let Some(e) = samples[0].ele {
        series.push((0.0, e));
    }
    for i in 1..samples.len() {
        dist_m += haversine_m(
            samples[i - 1].lat,
            samples[i - 1].lon,
            samples[i].lat,
            samples[i].lon,
        );
        if let Some(e) = samples[i].ele {
            series.push((dist_m, e));
        }
    }

    if series.len() < 2 {
        bail!("Not enough elevation data (need >= 2 points with elevation).");
    }

    let (min_d, max_d, min_e, max_e) = bounds_series(&series);

    // --- Layout / gutters so labels don’t clip ---
    let font_size = 12.0;
    let label_gap = 6.0;

    let left_gutter = if draw_labels { 48.0 } else { 0.0 };
    let bottom_gutter = if draw_labels { 18.0 } else { 0.0 };

    let default_plot_left = x + pad + left_gutter;
    let default_plot_right = x + w - pad;

    let plot_top = y + pad;
    let plot_bottom = y + h - pad - bottom_gutter;

    // Override plot width to match route outline width (embedded case)
    let mut plot_left = default_plot_left;
    let mut plot_right = default_plot_right;

    if let Some((l, r)) = plot_left_right {
        let min_left = x + pad + left_gutter;
        let max_right = x + w - pad;

        let candidate_left = l.max(min_left);
        let candidate_right = r.min(max_right);

        if candidate_right > candidate_left + 1.0 {
            plot_left = candidate_left;
            plot_right = candidate_right;
        }
    }

    let inner_w = (plot_right - plot_left).max(1.0);
    let inner_h = (plot_bottom - plot_top).max(1.0);

    // Axes align to plot area
    let axis_x = plot_left;
    let axis_y_top = plot_top;
    let axis_y_bottom = plot_bottom;
    let axis_x_right = plot_right;

    // Map to pixels
    let span_d = (max_d - min_d).abs().max(1e-12);
    let span_e = (max_e - min_e).abs().max(1e-12);

    let pts: Vec<Pt> = series
        .iter()
        .map(|(d, e)| Pt {
            x: plot_left + ((d - min_d) / span_d) * inner_w,
            y: plot_top + (1.0 - (e - min_e) / span_e) * inner_h,
        })
        .collect();

    // Axes (optional)
    let axes = if draw_axes {
        format!(
            r#"<g class="axes" stroke="{c}" stroke-width="1" fill="none">
  <line x1="{x:.2}" y1="{yt:.2}" x2="{x:.2}" y2="{yb:.2}" />
  <line x1="{x:.2}" y1="{yb:.2}" x2="{xr:.2}" y2="{yb:.2}" />
</g>"#,
            c = escape_attr(label_color),
            x = axis_x,
            yt = axis_y_top,
            yb = axis_y_bottom,
            xr = axis_x_right,
        )
    } else {
        String::new()
    };

    // Labels (optional)
    let labels = if draw_labels {
        let max_label_x = axis_x - label_gap;
        let max_label_y = axis_y_top;

        let min_label_x = axis_x - label_gap;
        let min_label_y = axis_y_bottom;

        let x0_label_x = axis_x;
        let x0_label_y = axis_y_bottom + font_size + 2.0;

        let x1_label_x = axis_x_right;
        let x1_label_y = x0_label_y;

        format!(
            r#"<g class="labels" font-size="{fs}" font-family="system-ui, -apple-system, sans-serif" fill="{c}">
  <text x="{max_x:.2}" y="{max_y:.2}" text-anchor="end" dominant-baseline="hanging">{max_e:.0} m</text>
  <text x="{min_x:.2}" y="{min_y:.2}" text-anchor="end" dominant-baseline="alphabetic">{min_e:.0} m</text>
  <text x="{x0_x:.2}" y="{x0_y:.2}" text-anchor="start" dominant-baseline="hanging">0 km</text>
  <text x="{x1_x:.2}" y="{x1_y:.2}" text-anchor="end" dominant-baseline="hanging">{dist:.1} km</text>
</g>"#,
            fs = font_size,
            c = escape_attr(label_color),
            max_x = max_label_x,
            max_y = max_label_y,
            min_x = min_label_x,
            min_y = min_label_y,
            x0_x = x0_label_x,
            x0_y = x0_label_y,
            x1_x = x1_label_x,
            x1_y = x1_label_y,
            max_e = max_e,
            min_e = min_e,
            dist = max_d / 1000.0,
        )
    } else {
        String::new()
    };

    // Fill under curve (optional, correct geometry)
    let fill_path = if fill != "none" {
        let mut d = String::new();
        if let (Some(first), Some(last)) = (pts.first(), pts.last()) {
            d.push_str(&format!("M {:.2} {:.2}", first.x, axis_y_bottom));
            d.push_str(&format!(" L {:.2} {:.2}", first.x, first.y));
            for p in &pts[1..] {
                d.push_str(&format!(" L {:.2} {:.2}", p.x, p.y));
            }
            d.push_str(&format!(" L {:.2} {:.2} Z", last.x, axis_y_bottom));
        }
        format!(
            r#"<path d="{d}" fill="{fill}" stroke="none" />"#,
            d = d,
            fill = escape_attr(fill)
        )
    } else {
        String::new()
    };

    // Elevation line (on top)
    let line = format!(
        r#"<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{sw}"
  stroke-linejoin="round" stroke-linecap="round" />"#,
        pts = pts_to_polyline(&pts),
        stroke = escape_attr(stroke),
        sw = sw
    );

    Ok(format!(
        r#"<g class="elevation">
  {axes}
  {labels}
  {fill_path}
  {line}
</g>"#,
        axes = axes,
        labels = labels,
        fill_path = fill_path,
        line = line
    ))
}

// -------------------------
// Route helpers
// -------------------------

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

fn px_x_span(pts: &[Pt]) -> (f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    for p in pts {
        min_x = min_x.min(p.x);
        max_x = max_x.max(p.x);
    }
    (min_x, max_x)
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

fn svg_teardrop(p: Pt, size: f64, fill: &str, stroke: &str, stroke_width: f64) -> String {
    // Teardrop pointing down, with the "tip" at (p.x, p.y)
    let r = size;
    let cy = p.y - r; // center of circular part

    format!(
        r#"<path d="
M {x:.2} {tip_y:.2}
C {x1:.2} {y1:.2}, {x2:.2} {y1:.2}, {x2:.2} {cy:.2}
A {r:.2} {r:.2} 0 1 1 {x0:.2} {cy:.2}
C {x0:.2} {y1:.2}, {x1:.2} {y1:.2}, {x:.2} {tip_y:.2}
Z"
fill="{fill}" stroke="{stroke}" stroke-width="{sw}" />"#,
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

// -------------------------
// Elevation helpers
// -------------------------

fn bounds_series(series: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;

    for (x, y) in series {
        min_x = min_x.min(*x);
        max_x = max_x.max(*x);
        min_y = min_y.min(*y);
        max_y = max_y.max(*y);
    }
    (min_x, max_x, min_y, max_y)
}

fn pts_to_polyline(pts: &[Pt]) -> String {
    let mut s = String::new();
    for (i, p) in pts.iter().enumerate() {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&format!("{:.2},{:.2}", p.x, p.y));
    }
    s
}

// Haversine distance (meters)
fn haversine_m(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let r = 6_371_000.0_f64;
    let dlat = (lat2 - lat1).to_radians();
    let dlon = (lon2 - lon1).to_radians();

    let a = (dlat / 2.0).sin().powi(2)
        + lat1.to_radians().cos() * lat2.to_radians().cos() * (dlon / 2.0).sin().powi(2);

    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
    r * c
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

// ---------------------
// Common escaping
// ---------------------

fn escape_attr(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

