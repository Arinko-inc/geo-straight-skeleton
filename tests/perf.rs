//! Timing of large inputs. Run with `cargo test --release --test perf -- --nocapture` to see the
//! numbers; in release builds a 5,000-vertex polygon must take less than a second.

use geo_straight_skeleton::Skeleton;
use std::time::{Duration, Instant};

fn rng(seed: u64) -> impl FnMut() -> f64 {
    let mut s = seed | 1;
    move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A circle of radius 100 with a 2% radial jitter (many reflex vertices).
fn noisy_circle(n: usize, seed: u64) -> Vec<[f64; 2]> {
    let mut r = rng(seed);
    (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64 + 0.3 * (r() - 0.5)) / n as f64;
            let rad = 100.0 * (1.0 + 0.02 * (2.0 * r() - 1.0));
            [rad * a.cos(), rad * a.sin()]
        })
        .collect()
}

/// A road-like strip: 3 m wide, gently curving, `n` vertices.
fn road(n: usize, seed: u64) -> Vec<[f64; 2]> {
    let mut r = rng(seed);
    let segs = n / 2 - 1;
    let mut c = vec![[0.0f64, 0.0f64]];
    let (mut ang, mut dang) = (0.0f64, 0.0f64);
    for _ in 0..segs {
        dang = 0.7 * dang + 0.03 * (r() - 0.5);
        ang += dang;
        let p = c[c.len() - 1];
        c.push([p[0] + 2.0 * ang.cos(), p[1] + 2.0 * ang.sin()]);
    }
    let mut left = Vec::new();
    let mut right = Vec::new();
    for i in 0..c.len() {
        let a = c[i.saturating_sub(1)];
        let b = c[(i + 1).min(c.len() - 1)];
        let d = [b[0] - a[0], b[1] - a[1]];
        let l = d[0].hypot(d[1]);
        let nrm = [-d[1] / l * 1.5, d[0] / l * 1.5];
        left.push([c[i][0] + nrm[0], c[i][1] + nrm[1]]);
        right.push([c[i][0] - nrm[0], c[i][1] - nrm[1]]);
    }
    left.reverse();
    right.extend(left);
    right
}

fn time(rings: &[Vec<[f64; 2]>], reps: usize) -> Duration {
    let mut best = Duration::MAX;
    for _ in 0..reps {
        let t = Instant::now();
        let sk = Skeleton::build(rings).expect("valid input");
        best = best.min(t.elapsed());
        assert_eq!(sk.faces.len(), rings.iter().map(|r| r.len()).sum::<usize>());
    }
    best
}

#[test]
fn large_inputs() {
    let release = !cfg!(debug_assertions);
    let reps = if release { 5 } else { 1 };
    let cases: Vec<(&str, Vec<Vec<[f64; 2]>>)> = vec![
        ("noisy circle, 1,000 vertices", vec![noisy_circle(1000, 1)]),
        ("noisy circle, 5,000 vertices", vec![noisy_circle(5000, 2)]),
        ("road strip, 2,000 vertices", vec![road(2000, 3)]),
        ("road strip, 5,000 vertices", vec![road(5000, 4)]),
    ];
    for (name, rings) in &cases {
        let t = time(rings, reps);
        eprintln!("{name}: {:.1} ms", t.as_secs_f64() * 1e3);
        if release {
            assert!(t < Duration::from_secs(1), "{name} took {t:?}");
        }
    }
}
