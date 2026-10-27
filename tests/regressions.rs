//! Regressions found by review (2026-10-27). Each input once produced a wrong `Ok`.

mod common;

use common::area as ring_area;
use geo_straight_skeleton::{Skeleton, SkeletonError};

/// A very sharp but positive input angle (1e-11 rad) must not be taken for a zero-width overlap:
/// the triangle has an incentre at about (1e12 - 5, 5) with time 5. If it cannot be computed,
/// the result must be `Numerical`, never a skeleton without nodes.
#[test]
fn sharp_input_angle_is_not_zero_width() {
    let tri = vec![[0.0, 0.0], [1e12, 0.0], [1e12, 10.0]];
    match Skeleton::build(&[tri]) {
        Ok(sk) => {
            let nodes: Vec<[f64; 3]> = sk.vertices.iter().copied().filter(|v| v[2] > 0.0).collect();
            assert_eq!(nodes.len(), 1, "{:?}", sk.vertices);
            let n = nodes[0];
            assert!((n[2] - 5.0).abs() < 1e-3 && (n[1] - 5.0).abs() < 1e-3 && (n[0] - (1e12 - 5.0)).abs() < 1e-2, "{n:?}");
            assert_eq!(sk.arcs().len(), 3);
        }
        Err(e) => assert_eq!(e, SkeletonError::Numerical),
    }
}

/// Points that input processing treats as duplicates (closer than 1e-12 of the extent) must still
/// join up in `offset`.
#[test]
fn offset_with_near_duplicate_points() {
    let e = 1e-12;
    let outer = vec![[0.0, 0.0], [e, 0.0], [4.0, 0.0], [4.0, e], [4.0, 4.0], [4.0 - e, 4.0], [0.0, 4.0], [0.0, 4.0 - e]];
    let sk = Skeleton::build(&[outer]).unwrap();
    let rings = sk.offset(0.5).unwrap();
    assert_eq!(rings.len(), 1, "{rings:?}");
    assert!((ring_area(&rings[0]) - 9.0).abs() < 1e-9, "{rings:?}");
}

/// The same for rings that touch at a vertex given as +0.0 in one ring and -0.0 in the other.
#[test]
fn offset_with_touching_rings_and_signed_zero() {
    let cases = [
        (vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]], vec![[-0.0, 0.0], [2.0, 5.0], [5.0, 2.0]]),
        (vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 5.0]], vec![[-0.0, 5.0], [3.0, 7.0], [3.0, 3.0]]),
        (vec![[0.0, -5.0], [10.0, -5.0], [10.0, 5.0], [0.0, 5.0], [-0.0, 0.0]], vec![[0.0, 0.0], [3.0, 2.0], [3.0, -2.0]]),
    ];
    for (outer, hole) in cases {
    let sk = Skeleton::build(&[outer, hole]).unwrap();
    for d in [0.1, 0.5, 1.0] {
        let rings = sk.offset(d).unwrap();
        let area: f64 = rings.iter().map(|r| ring_area(r)).sum();
        // the area of the faces above d, from the faces themselves
        let mut want = 0.0;
        for f in &sk.faces {
            let vs: Vec<[f64; 3]> = f.vertices.iter().map(|&i| sk.vertices[i]).collect();
            let mut clip: Vec<[f64; 2]> = Vec::new();
            for i in 0..vs.len() {
                let (p, q) = (vs[i], vs[(i + 1) % vs.len()]);
                if p[2] >= d {
                    clip.push([p[0], p[1]]);
                }
                if (p[2] >= d) != (q[2] >= d) {
                    let t = (d - p[2]) / (q[2] - p[2]);
                    clip.push([p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t]);
                }
            }
            if clip.len() >= 3 {
                want += ring_area(&clip);
            }
        }
        assert!((area - want).abs() < 1e-9, "offset {d}: {area} vs {want}: {rings:?}");
        // the hole has come off the outer ring: one ring goes round both
        assert_eq!(rings.len(), 1, "offset {d}: {rings:?}");
    }
    }
}

/// A small hole that grows under the offset must survive: offsets drop only collapsed contours,
/// not small ones.
#[test]
fn offset_keeps_small_hole() {
    let outer = vec![[0.0, 0.0], [1000.0, 0.0], [1000.0, 1000.0], [0.0, 1000.0]];
    let hole = vec![[500.0, 500.0], [500.0, 500.01], [500.01, 500.01], [500.01, 500.0]];
    let sk = Skeleton::build(&[outer, hole]).unwrap();
    let rings = sk.offset(0.001).unwrap();
    assert_eq!(rings.len(), 2, "{rings:?}");
    let mut areas: Vec<f64> = rings.iter().map(|r| ring_area(r)).collect();
    areas.sort_by(f64::total_cmp);
    assert!((areas[0] + 0.012 * 0.012).abs() < 1e-9, "{areas:?}");
    assert!((areas[1] - 999.998 * 999.998).abs() < 1e-6, "{areas:?}");
}
