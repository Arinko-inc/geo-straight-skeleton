//! Fixed expectations: the test table of the first user's specification (Algo-main
//! `docs/spec/04-water-skeleton.md` §4) plus error handling.

use geo_straight_skeleton::{Skeleton, SkeletonError};

fn nodes(sk: &Skeleton) -> Vec<[f64; 3]> {
    sk.vertices.iter().copied().filter(|v| v[2] > 0.0).collect()
}

fn has_node(sk: &Skeleton, x: f64, y: f64, t: f64) -> bool {
    nodes(sk).iter().any(|v| (v[0] - x).abs() < 1e-9 && (v[1] - y).abs() < 1e-9 && (v[2] - t).abs() < 1e-9)
}

fn ring(pts: &[(f64, f64)]) -> Vec<[f64; 2]> {
    pts.iter().map(|&(x, y)| [x, y]).collect()
}

/// V, F, number of interior arcs.
fn counts(sk: &Skeleton) -> (usize, usize, usize) {
    (sk.vertices.len(), sk.faces.len(), sk.arcs().len())
}

#[test]
fn rectangle() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (4., 0.), (4., 2.), (0., 2.)])]).unwrap();
    assert_eq!(counts(&sk), (6, 4, 5));
    assert!(has_node(&sk, 1., 1., 1.));
    assert!(has_node(&sk, 3., 1., 1.));
    // faces follow the input edges and start at the edge's end point
    assert_eq!(sk.faces[0].edge, (0, 1));
    assert_eq!(sk.faces[0].vertices[0], 1);
    assert_eq!(*sk.faces[0].vertices.last().unwrap(), 0);
}

#[test]
fn square_meets_in_one_point() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (2., 0.), (2., 2.), (0., 2.)])]).unwrap();
    assert_eq!(counts(&sk), (5, 4, 4));
    assert!(has_node(&sk, 1., 1., 1.));
}

#[test]
fn square_100m_has_no_ridge() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (100., 0.), (100., 100.), (0., 100.)])]).unwrap();
    assert_eq!(nodes(&sk).len(), 1);
    assert!(has_node(&sk, 50., 50., 50.));
}

#[test]
fn l_shape() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (4., 0.), (4., 2.), (2., 2.), (2., 4.), (0., 4.)])]).unwrap();
    assert_eq!(counts(&sk), (9, 6, 8));
    for (x, y) in [(1., 1.), (3., 1.), (1., 3.)] {
        assert!(has_node(&sk, x, y, 1.), "node ({x}, {y})");
    }
}

#[test]
fn t_shape() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (6., 0.), (6., 2.), (4., 2.), (4., 6.), (2., 6.), (2., 2.), (0., 2.)])])
        .unwrap();
    assert_eq!(counts(&sk), (12, 8, 11));
    for (x, y) in [(1., 1.), (5., 1.), (3., 1.), (3., 5.)] {
        assert!(has_node(&sk, x, y, 1.), "node ({x}, {y})");
    }
}

#[test]
fn u_shape() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (6., 0.), (6., 5.), (4., 5.), (4., 2.), (2., 2.), (2., 5.), (0., 5.)])])
        .unwrap();
    assert_eq!(counts(&sk), (12, 8, 11));
    for (x, y) in [(1., 1.), (5., 1.), (1., 4.), (5., 4.)] {
        assert!(has_node(&sk, x, y, 1.), "node ({x}, {y})");
    }
}

#[test]
fn square_with_hole() {
    let outer = ring(&[(0., 0.), (10., 0.), (10., 10.), (0., 10.)]);
    let hole = ring(&[(3., 3.), (3., 7.), (7., 7.), (7., 3.)]);
    let sk = Skeleton::build(&[outer.clone(), hole.clone()]).unwrap();
    assert_eq!(counts(&sk), (12, 8, 12));
    for (x, y) in [(1.5, 1.5), (1.5, 8.5), (8.5, 8.5), (8.5, 1.5)] {
        assert!(has_node(&sk, x, y, 1.5), "node ({x}, {y})");
    }
    // a counter-clockwise hole is flipped, not rejected
    let mut ccw = hole;
    ccw.reverse();
    let sk2 = Skeleton::build(&[outer, ccw]).unwrap();
    assert_eq!(counts(&sk2), (12, 8, 12));
}

#[test]
fn thin_channel() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (200., 0.), (200., 3.), (0., 3.)])]).unwrap();
    assert_eq!(counts(&sk), (6, 4, 5));
    assert!(has_node(&sk, 1.5, 1.5, 1.5));
    assert!(has_node(&sk, 198.5, 1.5, 1.5));
}

fn regular(n: usize, r: f64) -> Vec<[f64; 2]> {
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            [r * a.cos(), r * a.sin()]
        })
        .collect()
}

#[test]
fn regular_polygons_collapse_to_one_node() {
    for n in [3, 4, 5, 6, 8, 12, 16, 17, 32, 64, 128, 256] {
        let sk = Skeleton::build(&[regular(n, 10.0)]).unwrap();
        let nd = nodes(&sk);
        assert_eq!(nd.len(), 1, "regular {n}-gon: {} nodes", nd.len());
        let t = 10.0 * (std::f64::consts::PI / n as f64).cos();
        assert!((nd[0][2] - t).abs() < 1e-9 && nd[0][0].abs() < 1e-9 && nd[0][1].abs() < 1e-9, "{n}: {:?}", nd[0]);
        assert_eq!(sk.arcs().len(), n, "regular {n}-gon arcs");
    }
}

#[test]
fn triangle_incenter() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (4., 0.), (2., 3.)])]).unwrap();
    assert_eq!(counts(&sk), (4, 3, 3));
    let r = 12.0 / (2.0 * 13f64.sqrt() + 4.0);
    assert!(has_node(&sk, 2.0, r, r));
}

#[test]
fn orientation_closing_and_duplicates() {
    let base = Skeleton::build(&[ring(&[(0., 0.), (4., 0.), (4., 2.), (0., 2.)])]).unwrap();
    let variants = [
        ring(&[(0., 2.), (4., 2.), (4., 0.), (0., 0.)]),
        ring(&[(0., 0.), (4., 0.), (4., 2.), (0., 2.), (0., 0.)]),
        ring(&[(0., 0.), (4., 0.), (4., 0.), (4., 2.), (0., 2.)]),
        ring(&[(0., 2.), (4., 2.), (4., 2.), (4., 0.), (0., 0.), (0., 2.)]),
    ];
    for v in variants {
        let sk = Skeleton::build(std::slice::from_ref(&v)).unwrap();
        let n_in = v.len() - usize::from(v.first() == v.last());
        assert_eq!(sk.faces.len(), n_in, "{v:?}");
        let a = nodes(&sk);
        let b = nodes(&base);
        assert_eq!(a.len(), b.len(), "{v:?}");
        for p in &a {
            assert!(has_node(&base, p[0], p[1], p[2]), "{v:?}: {p:?}");
        }
        let area: f64 = sk.faces.iter().map(|f| face_area(&sk, &f.vertices)).sum();
        assert!((area - 8.0).abs() < 1e-12, "{v:?}: area {area}");
    }
}

fn face_area(sk: &Skeleton, vs: &[usize]) -> f64 {
    let n = vs.len();
    (0..n)
        .map(|i| {
            let p = sk.vertices[vs[i]];
            let q = sk.vertices[vs[(i + 1) % n]];
            p[0] * q[1] - q[0] * p[1]
        })
        .sum::<f64>()
        * 0.5
}

#[test]
fn collinear_vertex_gets_its_own_face() {
    let sk = Skeleton::build(&[ring(&[(0., 0.), (2., 0.), (4., 0.), (4., 2.), (0., 2.)])]).unwrap();
    assert_eq!(counts(&sk), (8, 5, 7));
    assert!(has_node(&sk, 2., 1., 1.));
    assert!(sk.arcs().iter().any(|&(a, b)| a == 1 && sk.vertices[b][0] == 2.0 && sk.vertices[b][1] == 1.0));
}

#[test]
fn invalid_inputs_are_errors() {
    let bowtie = ring(&[(0., 0.), (4., 4.), (4., 0.), (0., 4.)]);
    assert!(matches!(Skeleton::build(&[bowtie]), Err(SkeletonError::SelfIntersection { .. })));
    assert_eq!(Skeleton::build(&[ring(&[(0., 0.), (1., 1.)])]), Err(SkeletonError::TooFewVertices { ring: 0 }));
    assert_eq!(
        Skeleton::build(&[ring(&[(0., 0.), (1., 1.), (2., 2.), (3., 3.)])]),
        Err(SkeletonError::ZeroArea { ring: 0 })
    );
    assert_eq!(Skeleton::build(&[]), Err(SkeletonError::EmptyInput));
    assert_eq!(
        Skeleton::build(&[ring(&[(0., 0.), (1., f64::NAN), (0., 1.)])]),
        Err(SkeletonError::NonFinite { ring: 0, index: 1 })
    );
    assert_eq!(
        Skeleton::build(&[ring(&[(0., 0.), (1., 0.), (f64::INFINITY, 1.)])]),
        Err(SkeletonError::NonFinite { ring: 0, index: 2 })
    );
    // all points equal
    assert_eq!(Skeleton::build(&[ring(&[(1., 1.), (1., 1.), (1., 1.)])]), Err(SkeletonError::TooFewVertices { ring: 0 }));
    // a spike that folds back on its own edge
    let fold = ring(&[(0., 0.), (4., 0.), (2., 0.), (2., 2.)]);
    assert!(matches!(Skeleton::build(&[fold]), Err(SkeletonError::SelfIntersection { .. })));
    // a vertex touching another edge
    let touch = ring(&[(0., 0.), (4., 0.), (4., 4.), (2., 0.), (0., 4.)]);
    assert!(matches!(Skeleton::build(&[touch]), Err(SkeletonError::SelfIntersection { .. })));
    let outer = ring(&[(0., 0.), (10., 0.), (10., 10.), (0., 10.)]);
    let outside = ring(&[(20., 20.), (20., 21.), (21., 21.)]);
    assert_eq!(Skeleton::build(&[outer.clone(), outside]), Err(SkeletonError::HoleOutside { ring: 1 }));
    let crossing = ring(&[(5., 5.), (5., 15.), (6., 15.), (6., 5.)]);
    assert!(matches!(Skeleton::build(&[outer.clone(), crossing]), Err(SkeletonError::SelfIntersection { .. })));
    let big = ring(&[(1., 1.), (1., 9.), (9., 9.), (9., 1.)]);
    let small = ring(&[(3., 3.), (3., 4.), (4., 4.), (4., 3.)]);
    assert_eq!(Skeleton::build(&[outer, big, small]), Err(SkeletonError::NestedHole { ring: 2, outer: 1 }));
}

fn ellipse(n: usize, seed: u64) -> Vec<[f64; 2]> {
    // the spec's case 13: an ellipse with a 2% jitter, about 100 m across, in metres
    let mut s = seed;
    let mut rnd = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            let r = 1.0 + 0.02 * (2.0 * rnd() - 1.0);
            [120.0 * r * a.cos() + 512_345.0, 80.0 * r * a.sin() + 3_840_210.0]
        })
        .collect()
}

#[test]
fn jittered_ellipses() {
    for (n, seed) in [(500, 7u64), (2000, 8)] {
        let r = ellipse(n, seed);
        let sk = Skeleton::build(&[r]).unwrap();
        assert_eq!(sk.faces.len(), n);
    }
}

#[test]
fn rings_touching_at_a_vertex() {
    // a triangular hole touching the outer ring at (5, 0)
    let outer = ring(&[(0., 0.), (5., 0.), (10., 0.), (10., 10.), (0., 10.)]);
    let hole = ring(&[(5., 0.), (3., 4.), (7., 4.)]);
    let sk = Skeleton::build(&[outer.clone(), hole]).unwrap();
    assert_eq!(sk.faces.len(), 8);
    let area: f64 = sk.faces.iter().map(|f| face_area(&sk, &f.vertices)).sum();
    assert!((area - 92.0).abs() < 1e-9, "area {area}");
    // two holes sharing a vertex, and three holes meeting in one point
    let two = [
        ring(&[(0., 0.), (20., 0.), (20., 10.), (0., 10.)]),
        ring(&[(2., 2.), (2., 8.), (10., 5.)]),
        ring(&[(10., 5.), (18., 8.), (18., 2.)]),
    ];
    assert_eq!(Skeleton::build(&two).unwrap().faces.len(), 10);
    let three = [
        ring(&[(0., 0.), (20., 0.), (20., 20.), (0., 20.)]),
        ring(&[(10., 10.), (4., 14.), (4., 6.)]),
        ring(&[(10., 10.), (16., 6.), (16., 14.)]),
        ring(&[(10., 10.), (6., 16.), (14., 16.)]),
    ];
    assert_eq!(Skeleton::build(&three).unwrap().faces.len(), 13);
    // overlapping wedges at a shared vertex are rejected
    let overlap = [
        ring(&[(0., 0.), (20., 0.), (20., 20.), (0., 20.)]),
        ring(&[(10., 10.), (4., 14.), (4., 6.)]),
        ring(&[(10., 10.), (5., 13.), (5., 7.)]),
    ];
    assert!(matches!(Skeleton::build(&overlap), Err(SkeletonError::SelfIntersection { .. })));
    // touching an edge (not a vertex) is still an error
    let on_edge = [outer, ring(&[(4., 0.), (3., 4.), (7., 4.)])];
    assert!(matches!(Skeleton::build(&on_edge), Err(SkeletonError::SelfIntersection { .. })));
}

fn ring_area(r: &[[f64; 2]]) -> f64 {
    let n = r.len();
    (0..n).map(|i| r[i][0] * r[(i + 1) % n][1] - r[(i + 1) % n][0] * r[i][1]).sum::<f64>() * 0.5
}

#[test]
fn offsets() {
    let rect = Skeleton::build(&[ring(&[(0., 0.), (4., 0.), (4., 2.), (0., 2.)])]).unwrap();
    let r = rect.offset(0.5).unwrap();
    assert_eq!(r.len(), 1);
    assert_eq!(r[0].len(), 4);
    assert!((ring_area(&r[0]) - 3.0).abs() < 1e-9);
    for p in &r[0] {
        assert!(((p[0] - 0.5).abs() < 1e-9 || (p[0] - 3.5).abs() < 1e-9) && ((p[1] - 0.5).abs() < 1e-9 || (p[1] - 1.5).abs() < 1e-9));
    }
    assert!(rect.offset(1.0).unwrap().is_empty());
    assert!(rect.offset(2.0).unwrap().is_empty());
    let back = rect.offset(0.0).unwrap();
    assert_eq!(back.len(), 1);
    assert!((ring_area(&back[0]) - 8.0).abs() < 1e-9);

    let outer = ring(&[(0., 0.), (10., 0.), (10., 10.), (0., 10.)]);
    let hole = ring(&[(3., 3.), (3., 7.), (7., 7.), (7., 3.)]);
    let sk = Skeleton::build(&[outer, hole]).unwrap();
    let r = sk.offset(1.0).unwrap();
    assert_eq!(r.len(), 2);
    let mut areas: Vec<f64> = r.iter().map(|x| ring_area(x)).collect();
    areas.sort_by(f64::total_cmp);
    assert!((areas[0] + 36.0).abs() < 1e-9 && (areas[1] - 64.0).abs() < 1e-9, "{areas:?}");

    // the U splits into two legs
    let u = Skeleton::build(&[ring(&[(0., 0.), (6., 0.), (6., 5.), (4., 5.), (4., 2.), (2., 2.), (2., 5.), (0., 5.)])]).unwrap();
    assert_eq!(u.offset(0.5).unwrap().len(), 1);
    // at 0.9 the bottom (height 2) still connects; legs are 2 wide
    assert_eq!(u.offset(0.9).unwrap().len(), 1);
}
