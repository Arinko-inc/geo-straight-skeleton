//! Shared helpers for the oracle and property tests.
#![allow(dead_code)]

use serde_json::Value;
use geo_straight_skeleton::Skeleton;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub type P2 = [f64; 2];

pub struct Case {
    pub name: String,
    pub group: String,
    pub rings: Vec<Vec<P2>>,
    pub cgal: Option<Cgal>,
    /// CGAL on a copy whose touching hole vertices were moved by `delta` (only for inputs whose
    /// rings touch at a vertex, which CGAL's wrapper rejects)
    pub cgal_perturbed: Option<(Cgal, f64)>,
    pub cgal_status: String,
    pub cgal_message: String,
    /// power-of-two extent the oracle used for float32 (rounding bound = s * 2^-24 per coordinate)
    pub frame_s: f64,
}

pub struct Cgal {
    pub vertices: Vec<[f64; 3]>,
    pub refined: Vec<bool>,
    pub polygons: Vec<Vec<usize>>,
}

pub fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("corpus")
}

fn f(v: &Value) -> f64 {
    v.as_f64().expect("number")
}

pub fn load_case(path: &Path) -> Case {
    let text = std::fs::read_to_string(path).expect("read corpus file");
    let v: Value = serde_json::from_str(&text).expect("parse corpus json");
    let rings = v["rings"]
        .as_array()
        .expect("rings")
        .iter()
        .map(|r| r.as_array().expect("ring").iter().map(|p| [f(&p[0]), f(&p[1])]).collect())
        .collect();
    let parse = |cg: &Value| -> Option<Cgal> {
        if cg["status"].as_str() != Some("ok") {
            return None;
        }
        Some(Cgal {
            vertices: cg["vertices"]
                .as_array()
                .expect("vertices")
                .iter()
                .map(|p| [f(&p[0]), f(&p[1]), f(&p[2])])
                .collect(),
            refined: cg["refined"].as_array().expect("refined").iter().map(|b| b.as_bool().expect("bool")).collect(),
            polygons: cg["polygons"]
                .as_array()
                .expect("polygons")
                .iter()
                .map(|p| p.as_array().expect("poly").iter().map(|i| i.as_u64().expect("index") as usize).collect())
                .collect(),
        })
    };
    let cg = &v["cgal"];
    let status = cg["status"].as_str().unwrap_or("missing").to_string();
    let message = cg["message"].as_str().unwrap_or("").to_string();
    let cgal = parse(cg);
    let cgal_perturbed = parse(&v["cgal_perturbed"]).map(|c| (c, v["cgal_perturbed"]["delta"].as_f64().expect("delta")));
    let frame_s = v["frame"]["s"].as_f64().unwrap_or(1.0);
    Case {
        name: v["name"].as_str().expect("name").to_string(),
        group: v["group"].as_str().unwrap_or("").to_string(),
        rings,
        cgal,
        cgal_perturbed,
        cgal_status: status,
        cgal_message: message,
        frame_s,
    }
}

pub fn load_corpus() -> Vec<Case> {
    // the corpus lives in the repository only, not in the published package
    let Ok(dir) = std::fs::read_dir(corpus_dir()) else {
        return Vec::new();
    };
    let mut paths: Vec<PathBuf> = dir
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().map(|x| x == "json").unwrap_or(false))
        .collect();
    paths.sort();
    if let Ok(only) = std::env::var("SK_CASE") {
        let names: Vec<&str> = only.split(',').collect();
        paths.retain(|p| names.iter().any(|n| p.file_stem().map(|s| s == *n).unwrap_or(false)));
    }
    paths.iter().map(|p| load_case(p)).collect()
}

/// Signed area, computed relative to the first point (large coordinates would lose digits).
pub fn area(poly: &[P2]) -> f64 {
    let n = poly.len();
    if n == 0 {
        return 0.0;
    }
    let o = poly[0];
    let mut a = 0.0;
    for i in 0..n {
        let p = [poly[i][0] - o[0], poly[i][1] - o[1]];
        let q = [poly[(i + 1) % n][0] - o[0], poly[(i + 1) % n][1] - o[1]];
        a += p[0] * q[1] - q[0] * p[1];
    }
    0.5 * a
}

/// Rings without closing point and consecutive duplicates (for areas and boundary distances).
pub fn clean_rings(rings: &[Vec<P2>]) -> Vec<Vec<P2>> {
    rings
        .iter()
        .map(|r| {
            let mut out: Vec<P2> = Vec::new();
            for &p in r {
                if out.last() != Some(&p) {
                    out.push(p);
                }
            }
            while out.len() > 1 && out[0] == out[out.len() - 1] {
                out.pop();
            }
            out
        })
        .collect()
}

pub fn polygon_area(rings: &[Vec<P2>]) -> f64 {
    let rs = clean_rings(rings);
    let mut a = area(&rs[0]).abs();
    for h in &rs[1..] {
        a -= area(h).abs();
    }
    a
}

pub fn bbox_diag(rings: &[Vec<P2>]) -> f64 {
    let (mut x0, mut y0, mut x1, mut y1) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for r in rings {
        for p in r {
            x0 = x0.min(p[0]);
            y0 = y0.min(p[1]);
            x1 = x1.max(p[0]);
            y1 = y1.max(p[1]);
        }
    }
    (x1 - x0).hypot(y1 - y0)
}

pub fn face_points(sk: &Skeleton, k: usize) -> Vec<P2> {
    sk.faces[k].vertices.iter().map(|&i| [sk.vertices[i][0], sk.vertices[i][1]]).collect()
}

fn seg_intersect_proper(a: P2, b: P2, c: P2, d: P2, tol: f64) -> bool {
    let o = |p: P2, q: P2, r: P2| (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
    let la = (b[0] - a[0]).hypot(b[1] - a[1]);
    let lc = (d[0] - c[0]).hypot(d[1] - c[1]);
    let d1 = o(c, d, a) / lc.max(1e-300);
    let d2 = o(c, d, b) / lc.max(1e-300);
    let d3 = o(a, b, c) / la.max(1e-300);
    let d4 = o(a, b, d) / la.max(1e-300);
    ((d1 > tol && d2 < -tol) || (d1 < -tol && d2 > tol)) && ((d3 > tol && d4 < -tol) || (d3 < -tol && d4 > tol))
}

/// A face is simple if no two non-adjacent edges cross (touching within `tol` is allowed,
/// since coincident skeleton nodes can make faces weakly simple).
pub fn face_is_simple(pts: &[P2], tol: f64) -> bool {
    let n = pts.len();
    if n < 3 {
        return true;
    }
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        for j in i + 2..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let c = pts[j];
            let d = pts[(j + 1) % n];
            if seg_intersect_proper(a, b, c, d, tol) {
                return false;
            }
        }
    }
    true
}

/// A straight skeleton face is monotone along its input edge: projected onto the edge direction,
/// its boundary goes once up and once down. `pts` starts with the edge's end point and ends with
/// its start point.
pub fn face_is_monotone(pts: &[P2], tol: f64) -> bool {
    let n = pts.len();
    if n < 4 {
        return true;
    }
    let (a, b) = (pts[n - 1], pts[0]);
    let d = [b[0] - a[0], b[1] - a[1]];
    let l = d[0].hypot(d[1]);
    if l == 0.0 {
        return true;
    }
    let s: Vec<f64> = pts.iter().map(|p| ((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / l).collect();
    // direction changes along the cycle, ignoring steps shorter than tol
    let mut dirs: Vec<i8> = Vec::new();
    for i in 0..n {
        let ds = s[(i + 1) % n] - s[i];
        if ds > tol {
            dirs.push(1);
        } else if ds < -tol {
            dirs.push(-1);
        }
    }
    let changes = (0..dirs.len()).filter(|&i| dirs[i] != dirs[(i + 1) % dirs.len()]).count();
    changes <= 2
}

pub fn point_in_poly(p: P2, poly: &[P2]) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let a = poly[i];
        let b = poly[j];
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if p[0] < x {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

pub fn point_seg_dist(p: P2, a: P2, b: P2) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let l2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if l2 > 0.0 { (((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / l2).clamp(0.0, 1.0) } else { 0.0 };
    let q = [a[0] + ab[0] * t, a[1] + ab[1] * t];
    (p[0] - q[0]).hypot(p[1] - q[1])
}

/// Group points closer than `tol` (union-find over a hash grid). Returns a cluster id per point.
pub fn cluster(points: &[[f64; 3]], tol: f64) -> Vec<usize> {
    let n = points.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let cell = tol.max(1e-300) * 2.0;
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        let key = ((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64);
        grid.entry(key).or_default().push(i);
    }
    for (i, p) in points.iter().enumerate() {
        let (cx, cy) = ((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(list) = grid.get(&(cx + dx, cy + dy)) {
                    for &j in list {
                        if j <= i {
                            continue;
                        }
                        let q = points[j];
                        if (p[0] - q[0]).hypot(p[1] - q[1]) <= tol && (p[2] - q[2]).abs() <= tol {
                            let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                            if a != b {
                                parent[a.max(b)] = a.min(b);
                            }
                        }
                    }
                }
            }
        }
    }
    (0..n).map(|i| find(&mut parent, i)).collect()
}

/// Arcs between clusters: faces' consecutive vertex pairs with at least one end at time > 0,
/// mapped to clusters, without self-loops and duplicates.
pub fn cluster_arcs(vertices: &[[f64; 3]], faces: &[Vec<usize>], cl: &[usize]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for f in faces {
        let m = f.len();
        for i in 0..m {
            let a = f[i];
            let b = f[(i + 1) % m];
            if vertices[a][2] <= 0.0 && vertices[b][2] <= 0.0 {
                continue;
            }
            let (ca, cb) = (cl[a], cl[b]);
            if ca != cb {
                out.push((ca.min(cb), ca.max(cb)));
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// The structural properties every result must have. Returns a list of problems.
pub fn check_properties(rings: &[Vec<P2>], sk: &Skeleton) -> Vec<String> {
    let mut problems = Vec::new();
    let clean: Vec<Vec<P2>> = rings
        .iter()
        .map(|r| {
            let mut r = r.clone();
            if r.len() > 1 && r[0] == r[r.len() - 1] {
                r.pop();
            }
            r
        })
        .collect();
    let n_edges: usize = clean.iter().map(|r| r.len()).sum();
    if sk.faces.len() != n_edges {
        problems.push(format!("{} faces for {} edges", sk.faces.len(), n_edges));
        return problems;
    }
    let diag = bbox_diag(rings);
    let target = polygon_area(rings);
    let mut sum = 0.0;
    for k in 0..sk.faces.len() {
        let pts = face_points(sk, k);
        if pts.len() < 3 {
            continue;
        }
        let a = area(&pts);
        if a < -1e-9 * diag * diag {
            problems.push(format!("face {k} has negative area {a}"));
        }
        sum += a;
        if !face_is_simple(&pts, 1e-9 * diag) {
            problems.push(format!("face {k} is not simple"));
        }
        // 1e-8: strips of exactly constant width put every event within ~1e-11 of the same time,
        // and the nodes then carry rounding noise of a few 1e-9 of the extent
        if !face_is_monotone(&pts, 1e-7 * diag) {
            problems.push(format!("face {k} is not monotone along its edge"));
        }
    }
    if (sum - target).abs() > 1e-9 * target.max(diag * diag * 1e-6) {
        problems.push(format!("face areas sum to {sum}, polygon area is {target}"));
    }
    // time = distance to the supporting line of every incident face; time <= distance to boundary
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); sk.vertices.len()];
    for (k, f) in sk.faces.iter().enumerate() {
        for &v in &f.vertices {
            incident[v].push(k);
        }
    }
    let boundary: Vec<(P2, P2)> = clean_rings(rings)
        .iter()
        .flat_map(|r| (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()])))
        .collect();
    let tol = 1e-7 * diag;
    for (i, v) in sk.vertices.iter().enumerate() {
        if v[2] <= 0.0 {
            continue;
        }
        let p = [v[0], v[1]];
        for &k in &incident[i] {
            let (s, e) = sk.faces[k].edge;
            let a = [sk.vertices[s][0], sk.vertices[s][1]];
            let b = [sk.vertices[e][0], sk.vertices[e][1]];
            let d = [b[0] - a[0], b[1] - a[1]];
            let l = d[0].hypot(d[1]);
            if l == 0.0 {
                continue;
            }
            let dist = (d[0] * (p[1] - a[1]) - d[1] * (p[0] - a[0])) / l;
            if (dist - v[2]).abs() > tol {
                problems.push(format!("vertex {i}: time {} but distance {dist} to the line of face {k}", v[2]));
                break;
            }
        }
        if boundary.len() <= 20000 {
            let dmin = boundary.iter().map(|&(a, b)| point_seg_dist(p, a, b)).fold(f64::INFINITY, f64::min);
            if v[2] > dmin + tol {
                problems.push(format!("vertex {i}: time {} exceeds the distance {dmin} to the boundary", v[2]));
            }
        }
        if problems.len() > 20 {
            break;
        }
    }
    problems
}

/// Hash grid for "is there a point within tol" queries.
pub struct Grid {
    cell: f64,
    map: HashMap<(i64, i64), Vec<usize>>,
}

impl Grid {
    pub fn new(points: &[[f64; 3]], tol: f64) -> Grid {
        let cell = tol.max(1e-300) * 2.0;
        let mut map: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, p) in points.iter().enumerate() {
            map.entry(((p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64)).or_default().push(i);
        }
        Grid { cell, map }
    }

    /// Indices of all points within `tol` (in x and y only).
    pub fn all_within(&self, p: &[f64; 3], points: &[[f64; 3]], tol: f64) -> Vec<usize> {
        let (cx, cy) = ((p[0] / self.cell).floor() as i64, (p[1] / self.cell).floor() as i64);
        let reach = (tol / self.cell).ceil() as i64;
        let mut out = Vec::new();
        for dx in -reach..=reach {
            for dy in -reach..=reach {
                if let Some(list) = self.map.get(&(cx + dx, cy + dy)) {
                    for &j in list {
                        let q = points[j];
                        if (p[0] - q[0]).hypot(p[1] - q[1]) <= tol {
                            out.push(j);
                        }
                    }
                }
            }
        }
        out
    }

    pub fn any_within(&self, p: &[f64; 3], points: &[[f64; 3]], tol: f64) -> bool {
        let (cx, cy) = ((p[0] / self.cell).floor() as i64, (p[1] / self.cell).floor() as i64);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(list) = self.map.get(&(cx + dx, cy + dy)) {
                    for &j in list {
                        let q = points[j];
                        if (p[0] - q[0]).hypot(p[1] - q[1]) <= tol && (p[2] - q[2]).abs() <= tol {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }
}
