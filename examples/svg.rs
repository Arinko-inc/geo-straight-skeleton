//! Read a polygon from JSON and draw its straight skeleton as SVG.
//!
//! ```text
//! cargo run --example svg -- polygon.json > skeleton.svg
//! cargo run --example svg -- polygon.json --offsets 0.5,1,1.5 > skeleton.svg
//! cargo run --example svg -- tests/corpus/cross_long.json --cgal > cgal.svg
//! ```
//!
//! The JSON is either an array of rings `[[[x, y], ...], ...]` (ring 0 is the outer boundary, the
//! others are holes) or an object with such a `rings` field (the corpus files). With `--cgal`, the
//! CGAL result stored in a corpus file is drawn instead of this crate's.

use geo_straight_skeleton::Skeleton;
use serde_json::Value;
use std::fmt::Write as _;

fn rings_of(v: &Value) -> Vec<Vec<[f64; 2]>> {
    let r = if v.is_array() { v } else { &v["rings"] };
    r.as_array()
        .expect("an array of rings")
        .iter()
        .map(|ring| {
            ring.as_array()
                .expect("a ring")
                .iter()
                .map(|p| [p[0].as_f64().expect("x"), p[1].as_f64().expect("y")])
                .collect()
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(path) = args.iter().find(|a| !a.starts_with("--") && !a.contains(',') && a.parse::<f64>().is_err())
    else {
        eprintln!("usage: svg <polygon.json> [--offsets d1,d2,...] [--cgal]");
        std::process::exit(2);
    };
    let offsets: Vec<f64> = args
        .iter()
        .position(|a| a == "--offsets")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or_default();
    let json: Value = serde_json::from_str(&std::fs::read_to_string(path).expect("read the file")).expect("valid JSON");
    let rings = rings_of(&json);

    // what to draw: vertices [x, y, time] and faces (vertex index lists)
    let (vertices, faces, title): (Vec<[f64; 3]>, Vec<Vec<usize>>, String) = if args.iter().any(|a| a == "--cgal") {
        let cg = &json["cgal"];
        let vs = cg["vertices"].as_array().expect("no CGAL result in this file");
        (
            vs.iter().map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap(), p[2].as_f64().unwrap()]).collect(),
            cg["polygons"]
                .as_array()
                .unwrap()
                .iter()
                .map(|f| f.as_array().unwrap().iter().map(|i| i.as_u64().unwrap() as usize).collect())
                .collect(),
            "CGAL".to_string(),
        )
    } else {
        match Skeleton::build(&rings) {
            Ok(sk) => (sk.vertices.clone(), sk.faces.iter().map(|f| f.vertices.clone()).collect(), "geo-straight-skeleton".into()),
            Err(e) => {
                eprintln!("error: {e}");
                std::process::exit(1);
            }
        }
    };

    // frame
    let (mut x0, mut y0, mut x1, mut y1) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in rings.iter().flatten() {
        x0 = x0.min(p[0]);
        y0 = y0.min(p[1]);
        x1 = x1.max(p[0]);
        y1 = y1.max(p[1]);
    }
    let (w, margin) = (900.0, 20.0);
    let s = (w - 2.0 * margin) / (x1 - x0).max(y1 - y0).max(f64::MIN_POSITIVE);
    let h = (y1 - y0) * s + 2.0 * margin;
    let tx = |x: f64| margin + (x - x0) * s;
    let ty = |y: f64| h - margin - (y - y0) * s;
    let path_of = |pts: &mut dyn Iterator<Item = [f64; 2]>| {
        let mut d = String::new();
        for (i, p) in pts.enumerate() {
            let _ = write!(d, "{}{:.3},{:.3} ", if i == 0 { "M" } else { "L" }, tx(p[0]), ty(p[1]));
        }
        d.push('Z');
        d
    };

    let mut out = String::new();
    let _ = writeln!(out, r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h:.0}" viewBox="0 0 {w} {h:.3}">"#);
    let _ = writeln!(out, "<title>{title}</title>");
    let mut d = String::new();
    for r in &rings {
        d.push_str(&path_of(&mut r.iter().copied()));
    }
    let _ = writeln!(out, r##"<path d="{d}" fill="#eef4fb" fill-rule="evenodd" stroke="#1f4e8c" stroke-width="1.5"/>"##);
    // faces, shaded by the mean time of their vertices
    let tmax = vertices.iter().map(|v| v[2]).fold(0.0, f64::max).max(f64::MIN_POSITIVE);
    for f in &faces {
        if f.len() < 3 {
            continue;
        }
        let mean = f.iter().map(|&i| vertices[i][2]).sum::<f64>() / f.len() as f64;
        let light = 96.0 - 30.0 * mean / tmax;
        let d = path_of(&mut f.iter().map(|&i| [vertices[i][0], vertices[i][1]]));
        let _ = writeln!(out, r##"<path d="{d}" fill="hsl(212,60%,{light:.0}%)" stroke="#2b6cb0" stroke-width="0.6"/>"##);
    }
    if !offsets.is_empty() {
        if let Ok(sk) = Skeleton::build(&rings) {
            for &dist in &offsets {
                let mut d = String::new();
                for r in sk.offset(dist) {
                    d.push_str(&path_of(&mut r.into_iter()));
                }
                let _ = writeln!(out, r##"<path d="{d}" fill="none" stroke="#c05621" stroke-width="1" stroke-dasharray="4 3"/>"##);
            }
        }
    }
    for v in vertices.iter().filter(|v| v[2] > 0.0) {
        let _ = writeln!(out, r##"<circle cx="{:.3}" cy="{:.3}" r="1.6" fill="#1a365d"/>"##, tx(v[0]), ty(v[1]));
    }
    out.push_str("</svg>\n");
    print!("{out}");
}
