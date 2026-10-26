//! Comparison with CGAL's straight skeleton (npm `straight-skeleton` 2.0.1, see `oracle/gen.mjs`).
//!
//! For every corpus case CGAL could solve:
//! * per face (matched by its input edge): area within a relative 1e-6;
//! * the face areas sum to the polygon area;
//! * every skeleton vertex (time > 0) has a counterpart on the other side within
//!   1e-6 x bounding-box diagonal (x, y and time);
//! * the number of interior arcs matches, counted between merged vertices (vertices at the same
//!   position are merged on each side first: CGAL does not merge coincident nodes).
//!
//! Tolerances also allow for the float32 rounding of CGAL vertices that could not be refined
//! (`refined[i] == false` in the corpus, about `s * 2^-24` per coordinate).
//!
//! Inputs whose rings touch at a vertex are rejected by CGAL's wrapper (it returns null). For those
//! the corpus holds CGAL's result on a copy whose touching hole vertices were pulled into their
//! holes by `delta` (about 1e-6 of the extent); the comparison then allows `10 * delta` more.
//! Cases CGAL could not solve at all are checked by `tests/properties.rs` only.

mod common;
use common::*;
use skeleton_rs::Skeleton;

struct Report {
    name: String,
    problems: Vec<String>,
}

#[derive(Clone, Copy)]
enum Mode {
    Exact,
    Perturbed(f64),
}

fn perimeter(pts: &[P2]) -> f64 {
    (0..pts.len())
        .map(|i| {
            let a = pts[i];
            let b = pts[(i + 1) % pts.len()];
            (b[0] - a[0]).hypot(b[1] - a[1])
        })
        .sum()
}

fn compare(case: &Case, cg: &Cgal, mode: Mode) -> Vec<String> {
    let mut problems = Vec::new();
    let sk = match Skeleton::build(&case.rings) {
        Ok(s) => s,
        Err(e) => return vec![format!("build failed: {e}")],
    };
    let diag = bbox_diag(&case.rings);
    let f32_err = case.frame_s * 2f64.powi(-23);
    let tol = 1e-6 * diag;
    let extra = match mode {
        Mode::Exact => 0.0,
        Mode::Perturbed(d) => 10.0 * d,
    };
    let unrefined_any = cg.refined.iter().any(|r| !r);
    let vtol = tol + extra + if unrefined_any { f32_err } else { 0.0 };

    // --- faces by input edge: the CGAL face whose base edge has (nearly) the same endpoints
    let starts: Vec<[f64; 3]> = cg
        .polygons
        .iter()
        .map(|poly| {
            let s = cg.vertices[*poly.last().expect("face")];
            [s[0], s[1], 0.0]
        })
        .collect();
    let etol = 2.0 * extra + 1e-12 * diag;
    let grid = Grid::new(&starts, etol.max(1e-9 * diag));
    let poly_pts = |poly: &[usize], vs: &[[f64; 3]]| -> Vec<P2> { poly.iter().map(|&i| [vs[i][0], vs[i][1]]).collect() };
    let mut matched = 0;
    for (k, face) in sk.faces.iter().enumerate() {
        if face.vertices.len() < 3 {
            continue;
        }
        let s = sk.vertices[face.edge.0];
        let e = sk.vertices[face.edge.1];
        let found = grid.all_within(&[s[0], s[1], 0.0], &starts, etol.max(1e-12)).into_iter().find(|&ck| {
            let end = cg.vertices[cg.polygons[ck][0]];
            (end[0] - e[0]).hypot(end[1] - e[1]) <= etol
        });
        let Some(ck) = found else {
            problems.push(format!("face {k}: no CGAL face for edge {:?}", face.edge));
            continue;
        };
        matched += 1;
        let ours = area(&face_points(&sk, k));
        let cpoly = &cg.polygons[ck];
        let cpts = poly_pts(cpoly, &cg.vertices);
        let theirs = area(&cpts);
        let perim = perimeter(&cpts);
        let face_unrefined = cpoly.iter().any(|&i| !cg.refined[i]);
        let atol = 1e-6 * theirs.abs() + perim * (extra + if face_unrefined { f32_err } else { 1e-12 * diag });
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
    if matches!(mode, Mode::Exact) && (sum_cgal - target).abs() > 1e-6 * target {
        problems.push(format!("CGAL faces sum to {sum_cgal}, polygon area {target}"));
    }

    // --- skeleton vertices: each one needs a counterpart on the other side. With a perturbed
    // oracle, CGAL has tiny extra nodes next to the touching points; they match our input vertices.
    let inner = |vs: &[[f64; 3]]| -> Vec<[f64; 3]> { vs.iter().copied().filter(|v| v[2] > 0.0).collect() };
    let (ours_v, cg_v) = (inner(&sk.vertices), inner(&cg.vertices));
    let (ours_all, cg_all): (Vec<[f64; 3]>, Vec<[f64; 3]>) = match mode {
        Mode::Exact => (ours_v.clone(), cg_v.clone()),
        Mode::Perturbed(_) => (sk.vertices.clone(), cg.vertices.clone()),
    };
    let ours_grid = Grid::new(&ours_all, vtol);
    let cg_grid = Grid::new(&cg_all, vtol);
    let mut missing = 0;
    for v in &cg_v {
        if !ours_grid.any_within(v, &ours_all, vtol) {
            missing += 1;
            if missing <= 5 {
                problems.push(format!("CGAL vertex {:?} has no match in ours", v));
            }
        }
    }
    let mut unmatched = 0;
    for v in &ours_v {
        if !cg_grid.any_within(v, &cg_all, vtol) {
            unmatched += 1;
            if unmatched <= 5 {
                problems.push(format!("our vertex {:?} has no match in CGAL", v));
            }
        }
    }
    if missing + unmatched > 10 {
        problems.push(format!("{missing} CGAL vertices unmatched, {unmatched} of ours unmatched"));
    }

    // --- arcs between merged vertices
    let ours_cl = cluster(&sk.vertices, vtol);
    let cg_cl = cluster(&cg.vertices, vtol);

    // --- exact mode: every face visits the same merged vertices, in the same order, as CGAL's
    if matches!(mode, Mode::Exact) {
        let cg_grid_all = Grid::new(&cg.vertices, vtol);
        let to_cg = |i: usize| -> Option<usize> {
            let v = sk.vertices[i];
            cg_grid_all
                .all_within(&v, &cg.vertices, vtol)
                .into_iter()
                .find(|&j| (cg.vertices[j][2] - v[2]).abs() <= vtol)
                .map(|j| cg_cl[j])
        };
        let squash = |seq: Vec<usize>| -> Vec<usize> {
            let mut out: Vec<usize> = Vec::new();
            for c in seq {
                if out.last() != Some(&c) {
                    out.push(c);
                }
            }
            while out.len() > 1 && out.first() == out.last() {
                out.pop();
            }
            out
        };
        let mut bad = 0;
        for (k, face) in sk.faces.iter().enumerate() {
            if face.vertices.len() < 3 {
                continue;
            }
            let Some(ours_seq) = face.vertices.iter().map(|&i| to_cg(i)).collect::<Option<Vec<usize>>>() else {
                continue; // reported by the vertex check
            };
            let s = sk.vertices[face.edge.0];
            let e = sk.vertices[face.edge.1];
            let Some(ck) = grid.all_within(&[s[0], s[1], 0.0], &starts, etol.max(1e-12)).into_iter().find(|&ck| {
                let end = cg.vertices[cg.polygons[ck][0]];
                (end[0] - e[0]).hypot(end[1] - e[1]) <= etol
            }) else {
                continue;
            };
            let a = squash(ours_seq);
            let b = squash(cg.polygons[ck].iter().map(|&j| cg_cl[j]).collect());
            if a != b {
                bad += 1;
                if bad <= 3 {
                    problems.push(format!("face {k}: merged vertex sequence {a:?} vs CGAL {b:?}"));
                }
            }
        }
        if bad > 3 {
            problems.push(format!("{bad} faces differ from CGAL in their vertex sequence"));
        }
    }
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
    if corpus.is_empty() {
        eprintln!("oracle: no corpus in {} (it is not part of the published crate)", corpus_dir().display());
        return;
    }
    let known: Vec<String> = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("known-diffs.md"),
    )
    .map(|s| {
        s.lines()
            .filter_map(|l| l.strip_prefix("- case: "))
            .map(|l| l.split_whitespace().next().unwrap_or("").trim_matches('`').to_string())
            .collect()
    })
    .unwrap_or_default();
    let (mut exact, mut perturbed) = (0, 0);
    let mut reports: Vec<Report> = Vec::new();
    let mut skipped: Vec<String> = Vec::new();
    for case in &corpus {
        let problems = if let Some(cg) = &case.cgal {
            exact += 1;
            compare(case, cg, Mode::Exact)
        } else if let Some((cg, delta)) = &case.cgal_perturbed {
            perturbed += 1;
            compare(case, cg, Mode::Perturbed(*delta))
        } else {
            skipped.push(format!("{} (CGAL: {} {})", case.name, case.cgal_status, case.cgal_message));
            continue;
        };
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
        "oracle: {} cases; {exact} compared with CGAL, {perturbed} with CGAL on a perturbed copy; {} passed, {} differ ({} unexpected); {} not solvable by CGAL: {}",
        corpus.len(),
        exact + perturbed - reports.len(),
        reports.len(),
        unexpected,
        skipped.len(),
        skipped.join(", ")
    );
    assert_eq!(unexpected, 0, "cases differ from CGAL");
}
