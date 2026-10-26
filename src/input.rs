//! Input validation and preparation: closing points, duplicates, orientation, normalisation,
//! self-intersection and hole containment checks.

use crate::SkeletonError;
use crate::geom::*;

/// Points closer than this (in normalised units, the larger side of the bounding box being 1)
/// are treated as duplicates.
const DUP_EPS: f64 = 1e-12;
/// Rings with an area below this (normalised units) are rejected as degenerate.
const AREA_EPS: f64 = 1e-12;
/// Distance below which two segments are considered touching (normalised units).
const TOUCH_EPS: f64 = 1e-12;

/// One computed (non-degenerate) edge, mapped back to the input.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CompEdge {
    /// Global index of the input edge (= index of its first vertex in the flattened input).
    pub orig: usize,
    /// Output vertex index at the start of the edge (in the computed orientation).
    pub start: usize,
    /// Output vertex index at the end of the edge.
    pub end: usize,
}

pub(crate) struct Prepared {
    /// Cleaned, oriented rings in normalised coordinates. Ring 0 is counter-clockwise, holes clockwise.
    pub rings: Vec<Vec<V2>>,
    /// `edges[r][j]` is the edge from `rings[r][j]` to `rings[r][j + 1]`.
    pub edges: Vec<Vec<CompEdge>>,
    /// Input edges of zero length: (global edge index, start vertex, end vertex).
    pub zero_edges: Vec<(usize, usize, usize)>,
    /// All input points (closing points removed), in input order.
    pub points: Vec<[f64; 2]>,
    pub center: V2,
    pub scale: f64,
}

impl Prepared {
    pub(crate) fn to_world(&self, p: V2) -> [f64; 2] {
        [p[0] * self.scale + self.center[0], p[1] * self.scale + self.center[1]]
    }
}

struct RawRing {
    /// offset of the ring's first point in the flattened input
    offset: usize,
    len: usize,
}

pub(crate) fn prepare(rings: &[Vec<[f64; 2]>]) -> Result<Prepared, SkeletonError> {
    if rings.is_empty() {
        return Err(SkeletonError::EmptyInput);
    }
    // 1. closing points, finiteness, flattened point list
    let mut points: Vec<[f64; 2]> = Vec::new();
    let mut raw: Vec<RawRing> = Vec::with_capacity(rings.len());
    for (ri, ring) in rings.iter().enumerate() {
        let mut len = ring.len();
        if len > 1 && ring[0] == ring[len - 1] {
            len -= 1;
        }
        for (i, p) in ring[..len].iter().enumerate() {
            if !p[0].is_finite() || !p[1].is_finite() {
                return Err(SkeletonError::NonFinite { ring: ri, index: i });
            }
        }
        if len < 3 {
            return Err(SkeletonError::TooFewVertices { ring: ri });
        }
        raw.push(RawRing { offset: points.len(), len });
        points.extend_from_slice(&ring[..len]);
    }
    // 2. normalisation frame from the bounding box of everything
    let (mut x0, mut y0, mut x1, mut y1) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in &points {
        x0 = x0.min(p[0]);
        y0 = y0.min(p[1]);
        x1 = x1.max(p[0]);
        y1 = y1.max(p[1]);
    }
    let center = [0.5 * (x0 + x1), 0.5 * (y0 + y1)];
    let scale = (x1 - x0).max(y1 - y0);
    if !(scale > 0.0) || !scale.is_finite() {
        return Err(SkeletonError::ZeroArea { ring: 0 });
    }
    let inv = 1.0 / scale;
    let norm_pt = |p: [f64; 2]| -> V2 { [(p[0] - center[0]) * inv, (p[1] - center[1]) * inv] };

    // 3. duplicates (runs of consecutive equal points, cyclic), orientation
    let mut out_rings: Vec<Vec<V2>> = Vec::with_capacity(raw.len());
    let mut out_edges: Vec<Vec<CompEdge>> = Vec::with_capacity(raw.len());
    let mut zero_edges = Vec::new();
    let mut zero_area_rings: Vec<usize> = Vec::new();
    for (ri, rr) in raw.iter().enumerate() {
        let pts: Vec<V2> = (0..rr.len).map(|i| norm_pt(points[rr.offset + i])).collect();
        let m = rr.len;
        // find a start index that is not a duplicate of its predecessor
        let same = |a: usize, b: usize| dist2(pts[a], pts[b]) <= DUP_EPS * DUP_EPS;
        let start = (0..m).find(|&i| !same(i, (i + m - 1) % m));
        let Some(start) = start else {
            return Err(SkeletonError::TooFewVertices { ring: ri });
        };
        let mut runs: Vec<(usize, usize)> = Vec::new();
        for k in 0..m {
            let i = (start + k) % m;
            if k > 0 && same(i, (i + m - 1) % m) {
                runs.last_mut().expect("run exists").1 = i;
            } else {
                runs.push((i, i));
            }
        }
        if runs.len() < 3 {
            return Err(SkeletonError::TooFewVertices { ring: ri });
        }
        let cleaned: Vec<V2> = runs.iter().map(|&(a, _)| pts[a]).collect();
        let area2 = ring_area2(&cleaned);
        if area2.abs() * 0.5 <= AREA_EPS {
            // a figure-eight has zero signed area too: report it as a self-intersection
            if !all_collinear(&cleaned) {
                zero_area_rings.push(ri);
            } else {
                return Err(SkeletonError::ZeroArea { ring: ri });
            }
        }
        let want_ccw = ri == 0;
        let reversed = (area2 > 0.0) != want_ccw;
        // zero-length input edges inside runs
        for &(a, b) in &runs {
            let mut i = a;
            while i != b {
                let j = (i + 1) % m;
                let (s, e) = if reversed { (j, i) } else { (i, j) };
                zero_edges.push((rr.offset + i, rr.offset + s, rr.offset + e));
                i = j;
            }
        }
        let n = runs.len();
        let (ring, edges): (Vec<V2>, Vec<CompEdge>) = if !reversed {
            let edges = (0..n)
                .map(|j| {
                    let a = runs[j];
                    let b = runs[(j + 1) % n];
                    CompEdge { orig: rr.offset + a.1, start: rr.offset + a.1, end: rr.offset + b.0 }
                })
                .collect();
            (cleaned, edges)
        } else {
            // computed ring = runs in reverse order
            let rev: Vec<(usize, usize)> = runs.iter().rev().copied().collect();
            let ring: Vec<V2> = rev.iter().map(|&(a, _)| pts[a]).collect();
            let edges = (0..n)
                .map(|j| {
                    // computed edge from rev[j] to rev[j+1]; in input order rev[j+1] precedes rev[j]
                    let a = rev[j];
                    let b = rev[(j + 1) % n];
                    CompEdge { orig: rr.offset + b.1, start: rr.offset + a.0, end: rr.offset + b.1 }
                })
                .collect();
            (ring, edges)
        };
        out_rings.push(ring);
        out_edges.push(edges);
    }

    check_simple(&out_rings, &out_edges, &raw)?;
    if let Some(&ring) = zero_area_rings.first() {
        return Err(SkeletonError::ZeroArea { ring });
    }
    check_holes(&out_rings)?;

    Ok(Prepared { rings: out_rings, edges: out_edges, zero_edges, points, center, scale })
}

fn all_collinear(pts: &[V2]) -> bool {
    let a = pts[0];
    let far = pts.iter().copied().fold(a, |m, p| if dist2(p, a) > dist2(m, a) { p } else { m });
    let d = sub(far, a);
    let l = norm(d);
    if l == 0.0 {
        return true;
    }
    pts.iter().all(|&p| (cross(d, sub(p, a)) / l).abs() <= TOUCH_EPS)
}

#[inline]
fn orient(a: V2, b: V2, c: V2) -> f64 {
    cross(sub(b, a), sub(c, a))
}

fn point_seg_dist2(p: V2, a: V2, b: V2) -> f64 {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    let t = if l2 > 0.0 { (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    dist2(p, add(a, mul(ab, t)))
}

fn segs_touch(a: V2, b: V2, c: V2, d: V2) -> bool {
    let d1 = orient(c, d, a);
    let d2 = orient(c, d, b);
    let d3 = orient(a, b, c);
    let d4 = orient(a, b, d);
    if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0)) {
        return true;
    }
    let e2 = TOUCH_EPS * TOUCH_EPS;
    point_seg_dist2(a, c, d) <= e2
        || point_seg_dist2(b, c, d) <= e2
        || point_seg_dist2(c, a, b) <= e2
        || point_seg_dist2(d, a, b) <= e2
}

fn check_simple(rings: &[Vec<V2>], edges: &[Vec<CompEdge>], raw: &[RawRing]) -> Result<(), SkeletonError> {
    struct Seg {
        ring: usize,
        j: usize,
        a: V2,
        b: V2,
        xmin: f64,
        xmax: f64,
        ymin: f64,
        ymax: f64,
    }
    let mut segs: Vec<Seg> = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        let n = ring.len();
        for j in 0..n {
            let a = ring[j];
            let b = ring[(j + 1) % n];
            segs.push(Seg {
                ring: r,
                j,
                a,
                b,
                xmin: a[0].min(b[0]),
                xmax: a[0].max(b[0]),
                ymin: a[1].min(b[1]),
                ymax: a[1].max(b[1]),
            });
        }
    }
    let err = |s: &Seg, t: &Seg| {
        let ea = edges[s.ring][s.j].orig - raw[s.ring].offset;
        let eb = edges[t.ring][t.j].orig - raw[t.ring].offset;
        SkeletonError::SelfIntersection { ring_a: s.ring, edge_a: ea, ring_b: t.ring, edge_b: eb }
    };
    segs.sort_by(|p, q| p.xmin.total_cmp(&q.xmin));
    for i in 0..segs.len() {
        let s = &segs[i];
        for t in &segs[i + 1..] {
            if t.xmin > s.xmax + TOUCH_EPS {
                break;
            }
            if t.ymin > s.ymax + TOUCH_EPS || s.ymin > t.ymax + TOUCH_EPS {
                continue;
            }
            if s.ring == t.ring {
                let n = rings[s.ring].len();
                let next = |j: usize| (j + 1) % n;
                if next(s.j) == t.j || next(t.j) == s.j {
                    // adjacent: they share one endpoint; they must not fold back onto each other
                    let (first, second) = if next(s.j) == t.j { (s, t) } else { (t, s) };
                    let e2 = TOUCH_EPS * TOUCH_EPS;
                    // shared vertex is first.b == second.a
                    let folds = point_seg_dist2(second.b, first.a, first.b) <= e2
                        || point_seg_dist2(first.a, second.a, second.b) <= e2;
                    if folds {
                        return Err(err(s, t));
                    }
                    continue;
                }
            }
            if segs_touch(s.a, s.b, t.a, t.b) {
                return Err(err(s, t));
            }
        }
    }
    Ok(())
}

fn point_in_ring(p: V2, ring: &[V2]) -> bool {
    let n = ring.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let a = ring[i];
        let b = ring[j];
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

fn check_holes(rings: &[Vec<V2>]) -> Result<(), SkeletonError> {
    for h in 1..rings.len() {
        // no edges touch, so one vertex decides containment
        let p = rings[h][0];
        if !point_in_ring(p, &rings[0]) {
            return Err(SkeletonError::HoleOutside { ring: h });
        }
        for (g, other) in rings.iter().enumerate().skip(1) {
            if g != h && point_in_ring(p, other) {
                return Err(SkeletonError::NestedHole { ring: h, outer: g });
            }
        }
    }
    Ok(())
}
