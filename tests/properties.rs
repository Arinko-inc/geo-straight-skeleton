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

fn cases() -> u32 {
    std::env::var("PROPTEST_CASES").ok().and_then(|s| s.parse().ok()).unwrap_or(300)
}

/// Small deterministic generator for the seeded shapes below.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        ((self.next() * n as f64) as usize).min(n - 1)
    }
}

/// Boundary of a random simply connected set of grid cells (no pinches), counter-clockwise,
/// collinear points removed, plus square holes inside cells surrounded by other cells.
fn polyomino(seed: u64, cells: usize, hole_side: f64) -> Vec<Vec<P2>> {
    use std::collections::{HashMap, HashSet};
    let mut rng = Rng(seed | 1);
    let mut set: HashSet<(i64, i64)> = HashSet::new();
    let mut list = vec![(0i64, 0i64)];
    set.insert((0, 0));
    let dirs = [(1, 0), (-1, 0), (0, 1), (0, -1)];
    let ring8 = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];
    let mut guard = 0;
    while list.len() < cells && guard < 100_000 {
        guard += 1;
        let (x, y) = list[rng.below(list.len())];
        let (dx, dy) = dirs[rng.below(4)];
        let (nx, ny) = (x + dx, y + dy);
        if set.contains(&(nx, ny)) {
            continue;
        }
        let occ: Vec<bool> = ring8.iter().map(|&(a, b)| set.contains(&(nx + a, ny + b))).collect();
        let runs = (0..8).filter(|&i| occ[i] && !occ[(i + 7) % 8]).count();
        let pinch = (1..8).step_by(2).any(|i| occ[i] && !occ[i - 1] && !occ[(i + 1) % 8]);
        if runs != 1 || pinch {
            continue;
        }
        set.insert((nx, ny));
        list.push((nx, ny));
    }
    let mut next: HashMap<(i64, i64), (i64, i64)> = HashMap::new();
    for &(x, y) in &list {
        let e = [((x, y), (x + 1, y), (0, -1)), ((x + 1, y), (x + 1, y + 1), (1, 0)), ((x + 1, y + 1), (x, y + 1), (0, 1)), ((x, y + 1), (x, y), (-1, 0))];
        for (a, b, (ox, oy)) in e {
            if !set.contains(&(x + ox, y + oy)) {
                next.insert(a, b);
            }
        }
    }
    let start = *next.keys().min().expect("boundary");
    let mut ring = vec![start];
    let mut cur = next[&start];
    while cur != start {
        ring.push(cur);
        cur = next[&cur];
    }
    let n = ring.len();
    let outer: Vec<P2> = (0..n)
        .filter(|&i| {
            let p = ring[(i + n - 1) % n];
            let q = ring[i];
            let r = ring[(i + 1) % n];
            (q.0 - p.0) * (r.1 - q.1) - (q.1 - p.1) * (r.0 - q.0) != 0
        })
        .map(|i| [ring[i].0 as f64, ring[i].1 as f64])
        .collect();
    let mut rings = vec![outer];
    if hole_side > 0.0 {
        for &(x, y) in &list {
            if ring8.iter().all(|&(a, b)| set.contains(&(x + a, y + b))) && rng.next() < 0.5 {
                let (c, h) = ([x as f64 + 0.5, y as f64 + 0.5], 0.5 * hole_side);
                rings.push(vec![[c[0] - h, c[1] - h], [c[0] - h, c[1] + h], [c[0] + h, c[1] + h], [c[0] + h, c[1] - h]]);
            }
        }
    }
    rings
}

/// A road-like strip: a smooth random centreline buffered to `width` (mitred joins).
fn strip(seed: u64, segs: usize, seg_len: f64, width: f64, curv: f64) -> Vec<P2> {
    let mut rng = Rng(seed | 1);
    let mut pts = vec![[0.0, 0.0]];
    let (mut ang, mut dang) = (0.0f64, 0.0f64);
    for _ in 0..segs {
        dang = 0.7 * dang + curv * (rng.next() - 0.5);
        ang += dang;
        let p = pts[pts.len() - 1];
        pts.push([p[0] + seg_len * ang.cos(), p[1] + seg_len * ang.sin()]);
    }
    let h = 0.5 * width;
    let (mut left, mut right) = (Vec::new(), Vec::new());
    for i in 0..pts.len() {
        let dir = |a: P2, b: P2| {
            let d = [b[0] - a[0], b[1] - a[1]];
            let l = d[0].hypot(d[1]);
            [d[0] / l, d[1] / l]
        };
        let (n, m) = if i > 0 && i + 1 < pts.len() {
            let (d1, d2) = (dir(pts[i - 1], pts[i]), dir(pts[i], pts[i + 1]));
            let n1 = [-d1[1], d1[0]];
            let s = [n1[0] - d2[1], n1[1] + d2[0]];
            let l = s[0].hypot(s[1]);
            let n = [s[0] / l, s[1] / l];
            (n, 1.0 / (n[0] * n1[0] + n[1] * n1[1]))
        } else {
            let d = if i == 0 { dir(pts[0], pts[1]) } else { dir(pts[i - 1], pts[i]) };
            ([-d[1], d[0]], 1.0)
        };
        let p = pts[i];
        left.push([p[0] + n[0] * h * m, p[1] + n[1] * h * m]);
        right.push([p[0] - n[0] * h * m, p[1] - n[1] * h * m]);
    }
    left.reverse();
    right.extend(left);
    right
}

/// Insert collinear points on some edges.
fn with_collinear(ring: &[P2], seed: u64, every: usize) -> Vec<P2> {
    let mut rng = Rng(seed | 1);
    let mut out = Vec::new();
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        out.push(a);
        if every > 0 && i % every == 0 {
            let f = 0.2 + 0.6 * rng.next();
            out.push([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f]);
        }
    }
    out
}

fn transform(rings: &[Vec<P2>], scale: f64, dx: f64, dy: f64) -> Vec<Vec<P2>> {
    rings.iter().map(|r| r.iter().map(|p| [p[0] * scale + dx, p[1] * scale + dy]).collect()).collect()
}

/// Valid inputs only: some generated shapes are rejected (rounding can make them touch).
fn check_if_valid(rings: Vec<Vec<P2>>) -> Result<(), TestCaseError> {
    match Skeleton::build(&rings) {
        Err(skeleton_rs::SkeletonError::Numerical) => Err(TestCaseError::fail(format!("numerical failure\n{rings:?}"))),
        Err(_) => Ok(()),
        Ok(sk) => {
            let mut problems = check_properties(&rings, &sk);
            problems.extend(sample_cover(&rings, &sk, 99, 60));
            prop_assert!(problems.is_empty(), "{:?}\n{:?}", problems, rings);
            Ok(())
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: cases(), max_global_rejects: 1_000_000, failure_persistence: None, ..ProptestConfig::default() })]

    #[test]
    fn rectilinear_with_holes(seed in any::<u64>(), cells in 4usize..80, hole in prop::sample::select(vec![0.0, 0.5, 0.25, 0.8]),
                              scale in prop::sample::select(vec![1.0, 3.0, 0.01, 1000.0])) {
        let rings = transform(&polyomino(seed, cells, hole), scale, -30_000.0 * (scale > 100.0) as u8 as f64, 0.0);
        check_if_valid(rings)?;
    }

    #[test]
    fn road_strips(seed in any::<u64>(), segs in 2usize..120, width in prop::sample::select(vec![2.0, 3.0, 5.0]),
                   curv in prop::sample::select(vec![0.0, 0.02, 0.1, 0.3]), collinear in 0usize..4,
                   far in any::<bool>()) {
        let mut r = strip(seed, segs, 4.0, width, curv);
        if collinear > 0 { r = with_collinear(&r, seed ^ 7, collinear); }
        let rings = if far { transform(&[r], 1.0, -30_000.0, 120_000.0) } else { vec![r] };
        check_if_valid(rings)?;
    }

    #[test]
    fn holes_near_or_touching(
        angles in prop::collection::vec(0.0f64..std::f64::consts::TAU, 8..40),
        radii in prop::collection::vec(6.0f64..10.0, 1..40),
        k in 0usize..40,
        gap in prop::sample::select(vec![0.0, 1e-9, 1e-6, 1e-3, 0.1]),
        size in 0.05f64..0.9,
        spread in 0.1f64..0.9,
    ) {
        let outer = star_polygon(&angles, &radii);
        prop_assume!(outer.len() >= 8);
        // a triangular hole whose apex is at (or `gap` inside of) outer vertex k
        let n = outer.len();
        let v = outer[k % n];
        let a = outer[(k % n + n - 1) % n];
        let b = outer[(k % n + 1) % n];
        let ang = |p: P2| (p[1] - v[1]).atan2(p[0] - v[0]);
        let (aa, ab) = (ang(b), ang(a));
        let mut sector = ab - aa;
        while sector <= 0.0 { sector += std::f64::consts::TAU; }
        prop_assume!(sector < std::f64::consts::PI);
        let mid = aa + 0.5 * sector;
        let half = 0.5 * sector * spread;
        let len = size * ((a[0] - v[0]).hypot(a[1] - v[1])).min((b[0] - v[0]).hypot(b[1] - v[1]));
        let apex = [v[0] + gap * mid.cos(), v[1] + gap * mid.sin()];
        let p1 = [apex[0] + len * (mid - half).cos(), apex[1] + len * (mid - half).sin()];
        let p2 = [apex[0] + len * (mid + half).cos(), apex[1] + len * (mid + half).sin()];
        let hole = if gap == 0.0 { vec![v, p2, p1] } else { vec![apex, p2, p1] };
        check_if_valid(vec![outer, hole])?;
    }

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
