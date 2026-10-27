//! Inward offset polygons from the skeleton faces.

use crate::Skeleton;
use std::collections::HashMap;

/// A crossing of the level `d` on the face edge `a -> b` (vertex indices, canonical).
#[derive(Clone, Copy)]
struct Crossing {
    a: usize,
    b: usize,
    point: [f64; 2],
    /// position along the face's input edge
    along: f64,
    up: bool,
}

impl Skeleton {
    /// The polygon shrunk inward by `distance`: the cross-section of the skeleton's roof at height
    /// `distance`. Returns closed rings without a closing point, counter-clockwise for outer
    /// boundaries and clockwise for holes; several outer rings appear once the polygon splits.
    /// Rings may touch at a point when `distance` equals the time of a skeleton vertex.
    ///
    /// `distance <= 0` gives the input rings (oriented as above); a distance beyond the largest
    /// vertex time gives an empty list.
    ///
    /// ```
    /// use geo_straight_skeleton::Skeleton;
    ///
    /// let sk = Skeleton::build(&[vec![[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]]]).unwrap();
    /// let rings = sk.offset(0.5);
    /// assert_eq!(rings.len(), 1);
    /// assert_eq!(rings[0].len(), 4); // the rectangle (0.5, 0.5) - (3.5, 1.5)
    /// ```
    pub fn offset(&self, distance: f64) -> Vec<Vec<[f64; 2]>> {
        if !distance.is_finite() {
            return Vec::new();
        }
        // A level exactly at the time of a skeleton vertex (common with axis-parallel input) makes
        // rings touch there; cutting a relative 1e-12 lower separates them cleanly.
        let d = (distance * (1.0 - 1e-12)).max(f64::MIN_POSITIVE);
        // the same position can carry several input indices (rings touching, duplicates)
        let mut by_pos: HashMap<(u64, u64), usize> = HashMap::new();
        let canon: Vec<usize> = self
            .vertices
            .iter()
            .enumerate()
            .map(|(i, v)| {
                if v[2] > 0.0 {
                    i
                } else {
                    *by_pos.entry((v[0].to_bits(), v[1].to_bits())).or_insert(i)
                }
            })
            .collect();
        let above = |i: usize| self.vertices[i][2] >= d;
        // segments of the level set inside each face: from a crossing going down to one going up
        let mut segments: Vec<(Crossing, Crossing)> = Vec::new();
        for face in &self.faces {
            let vs = &face.vertices;
            let m = vs.len();
            if m < 3 {
                continue;
            }
            let (s, e) = (self.vertices[face.edge.0], self.vertices[face.edge.1]);
            let dir = [e[0] - s[0], e[1] - s[1]];
            let mut cs: Vec<Crossing> = Vec::new();
            for i in 0..m {
                let (u, w) = (vs[i], vs[(i + 1) % m]);
                let (au, aw) = (above(u), above(w));
                if au == aw {
                    continue;
                }
                // interpolate from the lower index, so the face on the other side of this arc
                // computes the identical point
                let (lo, hi) = if u < w { (u, w) } else { (w, u) };
                let (pl, ph) = (self.vertices[lo], self.vertices[hi]);
                let f = ((d - pl[2]) / (ph[2] - pl[2])).clamp(0.0, 1.0);
                let point = [pl[0] + (ph[0] - pl[0]) * f, pl[1] + (ph[1] - pl[1]) * f];
                let along = (point[0] - s[0]) * dir[0] + (point[1] - s[1]) * dir[1];
                cs.push(Crossing { a: canon[u], b: canon[w], point, along, up: aw });
            }
            // along the edge direction the level line leaves the superlevel set at "down"
            // crossings and enters it at "up" crossings, alternately
            cs.sort_by(|x, y| x.along.total_cmp(&y.along).then(x.up.cmp(&y.up)));
            let mut i = 0;
            while i + 1 < cs.len() {
                if !cs[i].up && cs[i + 1].up {
                    segments.push((cs[i], cs[i + 1]));
                    i += 2;
                } else {
                    i += 1;
                }
            }
        }
        // chain: a segment ending on face edge (a, b) continues in the face that has (b, a)
        let mut start_of: HashMap<(usize, usize), usize> = HashMap::with_capacity(segments.len());
        for (k, (from, _)) in segments.iter().enumerate() {
            start_of.insert((from.a, from.b), k);
        }
        let mut used = vec![false; segments.len()];
        let mut rings = Vec::new();
        // rings thinner than rounding (a cut at the very top of a ridge) are dropped
        let (mut x0, mut y0, mut x1, mut y1) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for v in &self.vertices {
            x0 = x0.min(v[0]);
            y0 = y0.min(v[1]);
            x1 = x1.max(v[0]);
            y1 = y1.max(v[1]);
        }
        let min_area2 = 2e-10 * ((x1 - x0).powi(2) + (y1 - y0).powi(2));
        for k0 in 0..segments.len() {
            if used[k0] {
                continue;
            }
            let mut ring: Vec<[f64; 2]> = Vec::new();
            let mut k = k0;
            loop {
                used[k] = true;
                let (from, to) = segments[k];
                if ring.last() != Some(&from.point) {
                    ring.push(from.point);
                }
                match start_of.get(&(to.b, to.a)) {
                    Some(&nk) if !used[nk] => k = nk,
                    _ => {
                        if ring.last() != Some(&to.point) && ring.first() != Some(&to.point) {
                            ring.push(to.point);
                        }
                        break;
                    }
                }
            }
            while ring.len() > 1 && ring.first() == ring.last() {
                ring.pop();
            }
            if ring.len() >= 3 && area2(&ring).abs() > min_area2 {
                rings.push(ring);
            }
        }
        rings
    }
}

fn area2(r: &[[f64; 2]]) -> f64 {
    let o = r[0];
    (0..r.len())
        .map(|i| {
            let p = r[i];
            let q = r[(i + 1) % r.len()];
            (p[0] - o[0]) * (q[1] - o[1]) - (q[0] - o[0]) * (p[1] - o[1])
        })
        .sum()
}
