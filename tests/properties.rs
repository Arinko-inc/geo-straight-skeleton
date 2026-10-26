//! Structural properties, checked on every corpus case (including those CGAL cannot solve) and on
//! random simple polygons:
//! * no panic, no error for valid input;
//! * one face per input edge;
//! * every face is simple and counter-clockwise, the faces add up to the polygon area and do not
//!   overlap (sampled points lie in exactly one face);
//! * every skeleton vertex lies at distance `time` from the supporting line of each face around it,
//!   and `time` never exceeds the Euclidean distance to the boundary.

mod common;
use common::*;
use proptest::prelude::*;
use skeleton_rs::Skeleton;

fn sample_cover(rings: &[Vec<P2>], sk: &Skeleton, seed: u64, samples: usize) -> Vec<String> {
    let rs = clean_rings(rings);
    let (mut x0, mut y0, mut x1, mut y1) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in &rs[0] {
        x0 = x0.min(p[0]);
        y0 = y0.min(p[1]);
        x1 = x1.max(p[0]);
        y1 = y1.max(p[1]);
    }
    let faces: Vec<Vec<P2>> = (0..sk.faces.len()).map(|k| face_points(sk, k)).filter(|f| f.len() >= 3).collect();
    let mut s = seed | 1;
    let mut rnd = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut problems = Vec::new();
    let mut tried = 0;
    let mut inside_count = 0;
    while inside_count < samples && tried < samples * 50 {
        tried += 1;
        let p = [x0 + (x1 - x0) * rnd(), y0 + (y1 - y0) * rnd()];
        let inside = point_in_poly(p, &rs[0]) && !rs[1..].iter().any(|h| point_in_poly(p, h));
        if !inside {
            continue;
        }
        inside_count += 1;
        let hits = faces.iter().filter(|f| point_in_poly(p, f)).count();
        if hits != 1 {
            problems.push(format!("point {p:?} lies in {hits} faces"));
            if problems.len() > 5 {
                break;
            }
        }
    }
    problems
}

#[test]
fn corpus_properties() {
    let corpus = load_corpus();
    let mut failures = Vec::new();
    for case in &corpus {
        let sk = match Skeleton::build(&case.rings) {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{}: build failed: {e}", case.name));
                continue;
            }
        };
        let mut problems = check_properties(&case.rings, &sk);
        let n_pts: usize = case.rings.iter().map(|r| r.len()).sum();
        let samples = if n_pts > 1000 { 100 } else { 300 };
        problems.extend(sample_cover(&case.rings, &sk, 0x9e37_79b9_7f4a_7c15 ^ n_pts as u64, samples));
        for p in problems.iter().take(6) {
            failures.push(format!("{} (CGAL: {}): {p}", case.name, case.cgal_status));
        }
    }
    for f in &failures {
        eprintln!("{f}");
    }
    assert!(failures.is_empty(), "{} property failures", failures.len());
}

/// Random star-shaped polygon: sorted angles with random radii (always simple).
fn star_polygon(angles: &[f64], radii: &[f64]) -> Vec<P2> {
    let mut a: Vec<f64> = angles.to_vec();
    a.sort_by(f64::total_cmp);
    a.dedup_by(|x, y| (*x - *y).abs() < 1e-6);
    a.iter().zip(radii.iter().cycle()).map(|(&t, &r)| [r * t.cos(), r * t.sin()]).collect()
}

fn check_random(rings: Vec<Vec<P2>>) -> Result<(), TestCaseError> {
    let sk = Skeleton::build(&rings).map_err(|e| TestCaseError::fail(format!("build failed: {e}")))?;
    let n: usize = rings.iter().map(|r| r.len()).sum();
    prop_assert_eq!(sk.faces.len(), n);
    let mut problems = check_properties(&rings, &sk);
    problems.extend(sample_cover(&rings, &sk, n as u64 * 7919, 60));
    prop_assert!(problems.is_empty(), "{:?}\n{:?}", problems, rings);
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 400, ..ProptestConfig::default() })]

    #[test]
    fn random_star_polygons(
        angles in prop::collection::vec(0.0f64..std::f64::consts::TAU, 3..60),
        radii in prop::collection::vec(1.0f64..10.0, 1..60),
    ) {
        let poly = star_polygon(&angles, &radii);
        prop_assume!(poly.len() >= 3);
        // the centre must be strictly inside: no angular gap of pi or more
        let mut a: Vec<f64> = poly.iter().map(|p| p[1].atan2(p[0])).collect();
        a.sort_by(f64::total_cmp);
        let mut gap = a[0] + std::f64::consts::TAU - a[a.len() - 1];
        for w in a.windows(2) { gap = gap.max(w[1] - w[0]); }
        prop_assume!(gap < std::f64::consts::PI - 1e-3);
        check_random(vec![poly])?;
    }

    #[test]
    fn random_grid_polygons(
        angles in prop::collection::vec(0usize..48, 3..40),
        radii in prop::collection::vec(2u32..12, 1..40),
    ) {
        // star-shaped polygons on an integer grid: many exactly parallel and collinear edges
        let mut a: Vec<usize> = angles.clone();
        a.sort();
        a.dedup();
        prop_assume!(a.len() >= 3);
        let poly: Vec<P2> = a.iter().zip(radii.iter().cycle()).map(|(&k, &r)| {
            let t = std::f64::consts::TAU * k as f64 / 48.0;
            [(r as f64 * t.cos()).round(), (r as f64 * t.sin()).round()]
        }).collect();
        let mut clean: Vec<P2> = Vec::new();
        for p in poly { if clean.last() != Some(&p) { clean.push(p); } }
        while clean.len() > 1 && clean[0] == clean[clean.len() - 1] { clean.pop(); }
        prop_assume!(clean.len() >= 3);
        // rounding can break simplicity or make the ring degenerate; only valid input is checked
        if let Ok(sk) = Skeleton::build(&[clean.clone()]) {
            let mut problems = check_properties(&[clean.clone()], &sk);
            problems.extend(sample_cover(&[clean.clone()], &sk, 17, 60));
            prop_assert!(problems.is_empty(), "{:?}\n{:?}", problems, clean);
        }
    }

    #[test]
    fn random_polygons_with_holes(
        angles in prop::collection::vec(0.0f64..std::f64::consts::TAU, 8..40),
        radii in prop::collection::vec(6.0f64..10.0, 1..40),
        hole_n in 3usize..10,
        hole_r in 1.0f64..4.0,
        hole_phase in 0.0f64..1.0,
    ) {
        let outer = star_polygon(&angles, &radii);
        prop_assume!(outer.len() >= 8);
        let mut a: Vec<f64> = outer.iter().map(|p| p[1].atan2(p[0])).collect();
        a.sort_by(f64::total_cmp);
        let mut gap = a[0] + std::f64::consts::TAU - a[a.len() - 1];
        for w in a.windows(2) { gap = gap.max(w[1] - w[0]); }
        prop_assume!(gap < 1.0);
        let hole: Vec<P2> = (0..hole_n).rev().map(|i| {
            let t = hole_phase + std::f64::consts::TAU * i as f64 / hole_n as f64;
            [hole_r * t.cos(), hole_r * t.sin()]
        }).collect();
        check_random(vec![outer, hole])?;
    }
}
