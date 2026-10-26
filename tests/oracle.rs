//! Comparison with CGAL's straight skeleton (npm `straight-skeleton` 2.0.1, see `oracle/gen.mjs`).
//!
//! For every corpus case CGAL could solve:
//! * per face (matched by its input edge): area within a relative 1e-6;
//! * the face areas sum to the polygon area;
//! * skeleton vertices (time > 0): after merging vertices at the same position on each side, every
//!   vertex has a counterpart on the other side within 1e-6 x bounding-box diagonal (x, y and time);
//! * the number of interior arcs matches (counted between merged vertices).
//!
//! Tolerances also allow for the float32 rounding of CGAL vertices that could not be refined
//! (`refined[i] == false` in the corpus, about `s * 2^-24` per coordinate).
//! Cases CGAL could not solve are checked by `tests/properties.rs`.

mod common;
use common::*;
use skeleton_rs::Skeleton;
use std::collections::HashMap;

struct Report {
    name: String,
    problems: Vec<String>,
}

fn compare(case: &Case) -> Vec<String> {
    let cg = case.cgal.as_ref().expect("cgal result");
    let mut problems = Vec::new();
    let sk = match Skeleton::build(&case.rings) {
        Ok(s) => s,
        Err(e) => return vec![format!("build failed: {e}")],
    };
    let diag = bbox_diag(&case.rings);
    let f32_err = case.frame_s * 2f64.powi(-23);
    let tol = 1e-6 * diag;
    let unrefined_any = cg.refined.iter().any(|r| !r);
    let vtol = tol + if unrefined_any { f32_err } else { 0.0 };

    // --- faces by input edge (positions of the two endpoints)
    let key = |p: [f64; 3]| (p[0].to_bits(), p[1].to_bits());
    let mut cg_faces: HashMap<((u64, u64), (u64, u64)), usize> = HashMap::new();
    for (k, poly) in cg.polygons.iter().enumerate() {
        let end = cg.vertices[poly[0]];
        let start = cg.vertices[*poly.last().expect("face")];
        cg_faces.insert((key(start), key(end)), k);
    }
    let poly_pts = |poly: &[usize], vs: &[[f64; 3]]| -> Vec<P2> { poly.iter().map(|&i| [vs[i][0], vs[i][1]]).collect() };
    let mut matched = 0;
    for (k, face) in sk.faces.iter().enumerate() {
        if face.vertices.len() < 3 {
            continue;
        }
        let s = sk.vertices[face.edge.0];
        let e = sk.vertices[face.edge.1];
        let Some(&ck) = cg_faces.get(&(key(s), key(e))) else {
            problems.push(format!("face {k}: no CGAL face for edge {:?}", face.edge));
            continue;
        };
        matched += 1;
        let ours = area(&face_points(&sk, k));
        let cpoly = &cg.polygons[ck];
        let theirs = area(&poly_pts(cpoly, &cg.vertices));
        let perim: f64 = {
            let pts = poly_pts(cpoly, &cg.vertices);
            (0..pts.len())
                .map(|i| {
                    let a = pts[i];
                    let b = pts[(i + 1) % pts.len()];
                    (b[0] - a[0]).hypot(b[1] - a[1])
                })
                .sum()
        };
        let face_unrefined = cpoly.iter().any(|&i| !cg.refined[i]);
        let atol = 1e-6 * theirs.abs() + perim * if face_unrefined { f32_err } else { 1e-12 * diag };
        if (ours - theirs).abs() > atol {
            problems.push(format!(
                "face {k} edge {:?}: area {ours} vs CGAL {theirs} (rel {:.2e})",
                face.edge,
                (ours - theirs).abs() / theirs.abs().max(1e-300)
            ));
        }
    }
    if matched != cg.polygons.len() {
        problems.push(format!("matched {matched} faces, CGAL has {}", cg.polygons.len()));
    }
    let target = polygon_area(&case.rings);
    let sum_ours: f64 = (0..sk.faces.len()).map(|k| area(&face_points(&sk, k))).sum();
    let sum_cgal: f64 = cg.polygons.iter().map(|p| area(&poly_pts(p, &cg.vertices))).sum();
    if (sum_ours - target).abs() > 1e-9 * target {
        problems.push(format!("our faces sum to {sum_ours}, polygon area {target}"));
    }
    if (sum_cgal - target).abs() > 1e-6 * target {
        problems.push(format!("CGAL faces sum to {sum_cgal}, polygon area {target}"));
    }

    // --- skeleton vertices: each one needs a counterpart on the other side
    let ours_cl = cluster(&sk.vertices, vtol);
    let cg_cl = cluster(&cg.vertices, vtol);
    let inner = |vs: &[[f64; 3]]| -> Vec<[f64; 3]> { vs.iter().copied().filter(|v| v[2] > 0.0).collect() };
    let ours_v = inner(&sk.vertices);
    let cg_v = inner(&cg.vertices);
    let ours_grid = Grid::new(&ours_v, vtol);
    let cg_grid = Grid::new(&cg_v, vtol);
    let near = |p: &[f64; 3], set: &[[f64; 3]]| {
        let grid = if std::ptr::eq(set, ours_v.as_slice()) { &ours_grid } else { &cg_grid };
        grid.any_within(p, set, vtol)
    };
    let mut missing = 0;
    for v in &cg_v {
        if !near(v, &ours_v) {
            missing += 1;
            if missing <= 5 {
                problems.push(format!("CGAL vertex {:?} has no match in ours", v));
            }
        }
    }
    let mut extra = 0;
    for v in &ours_v {
        if !near(v, &cg_v) {
            extra += 1;
            if extra <= 5 {
                problems.push(format!("our vertex {:?} has no match in CGAL", v));
            }
        }
    }
    if missing + extra > 10 {
        problems.push(format!("{missing} CGAL vertices unmatched, {extra} of ours unmatched"));
    }

    // --- arcs between merged vertices
    let ours_faces: Vec<Vec<usize>> = sk.faces.iter().map(|f| f.vertices.clone()).collect();
    let a_ours = cluster_arcs(&sk.vertices, &ours_faces, &ours_cl).len();
    let a_cg = cluster_arcs(&cg.vertices, &cg.polygons, &cg_cl).len();
    if a_ours != a_cg {
        problems.push(format!("{a_ours} arcs vs CGAL {a_cg} (raw arcs ours {})", sk.arcs().len()));
    }
    problems
}

#[test]
fn oracle_corpus() {
    let corpus = load_corpus();
    assert!(!corpus.is_empty(), "empty corpus");
    let known: Vec<String> = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("known-diffs.md"),
    )
    .map(|s| {
        s.lines()
            .filter_map(|l| l.strip_prefix("- case: "))
            .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
            .collect()
    })
    .unwrap_or_default();
    let mut compared = 0;
    let mut reports: Vec<Report> = Vec::new();
    let mut skipped: Vec<String> = Vec::new();
    for case in &corpus {
        if case.cgal.is_none() {
            skipped.push(format!("{} (CGAL: {} {})", case.name, case.cgal_status, case.cgal_message));
            continue;
        }
        compared += 1;
        let problems = compare(case);
        if !problems.is_empty() {
            reports.push(Report { name: case.name.clone(), problems });
        }
    }
    let mut unexpected = 0;
    for r in &reports {
        let is_known = known.contains(&r.name);
        if !is_known {
            unexpected += 1;
        }
        eprintln!("{} {}:", if is_known { "KNOWN" } else { "FAIL " }, r.name);
        for p in r.problems.iter().take(12) {
            eprintln!("    {p}");
        }
    }
    eprintln!(
        "oracle: {compared} compared, {} passed, {} differ ({} unexpected); {} not solvable by CGAL: {}",
        compared - reports.len(),
        reports.len(),
        unexpected,
        skipped.len(),
        skipped.join(", ")
    );
    assert_eq!(unexpected, 0, "cases differ from CGAL");
}
