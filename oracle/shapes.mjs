import { readdirSync, readFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

const inputsDir = join(dirname(fileURLToPath(import.meta.url)), 'inputs');

// Input generators for the oracle corpus. Every generator is deterministic (fixed seeds).
// A case is {name, group, note?, rings}. Rings carry no closing point unless the case says so;
// orientation is whatever the generator produces (gen.mjs fixes it before calling CGAL).

export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const TAU = Math.PI * 2;

function rect(x0, y0, w, h) {
  return [[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]];
}

function rotate(ring, ang, cx = 0, cy = 0) {
  const c = Math.cos(ang), s = Math.sin(ang);
  return ring.map(([x, y]) => [cx + (x - cx) * c - (y - cy) * s, cy + (x - cx) * s + (y - cy) * c]);
}

function translate(ring, dx, dy) {
  return ring.map(([x, y]) => [x + dx, y + dy]);
}

function regular(n, r, phase = 0, cx = 0, cy = 0) {
  const out = [];
  for (let i = 0; i < n; i++) {
    const a = phase + (TAU * i) / n;
    out.push([cx + r * Math.cos(a), cy + r * Math.sin(a)]);
  }
  return out;
}

function star(n, r1, r2, phase = 0) {
  const out = [];
  for (let i = 0; i < 2 * n; i++) {
    const a = phase + (Math.PI * i) / n;
    const r = i % 2 === 0 ? r1 : r2;
    out.push([r * Math.cos(a), r * Math.sin(a)]);
  }
  return out;
}

// Random star-shaped polygon: sorted random angles, random radii.
function randomStar(n, seed, rmin = 0.3, rmax = 1.0, scale = 10) {
  const rnd = mulberry32(seed);
  const angs = [];
  for (let i = 0; i < n; i++) angs.push(rnd() * TAU);
  angs.sort((a, b) => a - b);
  // keep angles distinct
  for (let i = 1; i < n; i++) if (angs[i] - angs[i - 1] < 1e-9) angs[i] = angs[i - 1] + 1e-9;
  return angs.map((a) => {
    const r = scale * (rmin + (rmax - rmin) * rnd());
    return [r * Math.cos(a), r * Math.sin(a)];
  });
}

// Noisy circle: n vertices, radius jitter `jit` (fraction), for the large cases.
function noisyCircle(n, seed, r = 100, jit = 0.02, aspect = 1) {
  const rnd = mulberry32(seed);
  const out = [];
  for (let i = 0; i < n; i++) {
    const a = (TAU * (i + 0.3 * (rnd() - 0.5))) / n;
    const rr = r * (1 + jit * (2 * rnd() - 1));
    out.push([rr * Math.cos(a) * aspect, rr * Math.sin(a)]);
  }
  return out;
}

function segIntersect(a, b, c, d) {
  const o = (p, q, r) => (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
  const d1 = o(c, d, a), d2 = o(c, d, b), d3 = o(a, b, c), d4 = o(a, b, d);
  return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

// Random simple polygon by 2-opt untangling of a random permutation of random points.
function random2opt(n, seed, scale = 10) {
  const rnd = mulberry32(seed);
  const pts = [];
  for (let i = 0; i < n; i++) pts.push([scale * rnd(), scale * rnd()]);
  let poly = pts.slice();
  for (let iter = 0; iter < 100000; iter++) {
    let changed = false;
    for (let i = 0; i < n && !changed; i++) {
      for (let j = i + 2; j < n && !changed; j++) {
        if (i === 0 && j === n - 1) continue;
        const a = poly[i], b = poly[i + 1], c = poly[j], d = poly[(j + 1) % n];
        if (segIntersect(a, b, c, d)) {
          // reverse i+1..j
          const mid = poly.slice(i + 1, j + 1).reverse();
          poly = [...poly.slice(0, i + 1), ...mid, ...poly.slice(j + 1)];
          changed = true;
        }
      }
    }
    if (!changed) break;
  }
  return poly;
}

// Random rectilinear polygon: boundary of a random simply-connected polyomino on a grid.
function randomPolyomino(cells, seed, unit = 1) {
  const rnd = mulberry32(seed);
  const key = (x, y) => `${x},${y}`;
  const set = new Set([key(0, 0)]);
  const list = [[0, 0]];
  const dirs = [[1, 0], [-1, 0], [0, 1], [0, -1]];
  // grow, rejecting cells that would create a hole (keep the complement connected approximately by
  // requiring the new cell to have exactly one occupied 4-neighbour run)
  let guard = 0;
  while (list.length < cells && guard++ < 100000) {
    const [x, y] = list[Math.floor(rnd() * list.length)];
    const [dx, dy] = dirs[Math.floor(rnd() * 4)];
    const nx = x + dx, ny = y + dy;
    if (set.has(key(nx, ny))) continue;
    // ring of 8 neighbours must have exactly one occupied run, to avoid holes and pinches
    const ring = [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]].map(([a, b]) => set.has(key(nx + a, ny + b)));
    let runs = 0;
    for (let i = 0; i < 8; i++) if (ring[i] && !ring[(i + 7) % 8]) runs++;
    // diagonal-only contacts create pinch points
    let bad = false;
    for (let i = 1; i < 8; i += 2) if (ring[i] && !ring[i - 1] && !ring[(i + 1) % 8]) bad = true;
    if (runs !== 1 || bad) continue;
    set.add(key(nx, ny));
    list.push([nx, ny]);
  }
  // boundary edges (directed CCW around cells), then chain
  const edges = new Map();
  for (const [x, y] of list) {
    const e = [[[x, y], [x + 1, y], 0, -1], [[x + 1, y], [x + 1, y + 1], 1, 0], [[x + 1, y + 1], [x, y + 1], 0, 1], [[x, y + 1], [x, y], -1, 0]];
    for (const [a, b, ox, oy] of e) {
      if (set.has(key(x + ox, y + oy))) continue;
      edges.set(key(a[0], a[1]), b);
    }
  }
  const start = edges.keys().next().value;
  const ring = [];
  let cur = start;
  let g = 0;
  do {
    const [cx, cy] = cur.split(',').map(Number);
    ring.push([cx, cy]);
    const nb = edges.get(cur);
    cur = key(nb[0], nb[1]);
  } while (cur !== start && g++ < 100000);
  // drop collinear vertices
  const out = [];
  for (let i = 0; i < ring.length; i++) {
    const p = ring[(i + ring.length - 1) % ring.length], q = ring[i], r = ring[(i + 1) % ring.length];
    const cr = (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]);
    if (cr !== 0) out.push([q[0] * unit, q[1] * unit]);
  }
  return out;
}

// A road-like strip: buffer a smooth random centerline by half-width w/2 (miter joins).
function road(nSeg, seed, segLen = 10, width = 3, curv = 0.08, opts = {}) {
  const rnd = mulberry32(seed);
  const pts = [[0, 0]];
  let ang = 0;
  let dang = 0;
  for (let i = 0; i < nSeg; i++) {
    dang = 0.7 * dang + curv * (rnd() - 0.5);
    ang += dang;
    const [x, y] = pts[pts.length - 1];
    pts.push([x + segLen * Math.cos(ang), y + segLen * Math.sin(ang)]);
  }
  const h = width / 2;
  const left = [], right = [];
  for (let i = 0; i < pts.length; i++) {
    const p = pts[i];
    const a = i > 0 ? pts[i - 1] : null, b = i < pts.length - 1 ? pts[i + 1] : null;
    const dir = (u, v) => { const dx = v[0] - u[0], dy = v[1] - u[1]; const l = Math.hypot(dx, dy); return [dx / l, dy / l]; };
    let n;
    let m = 1;
    if (a && b) {
      const d1 = dir(a, p), d2 = dir(p, b);
      const n1 = [-d1[1], d1[0]], n2 = [-d2[1], d2[0]];
      n = [n1[0] + n2[0], n1[1] + n2[1]];
      const l = Math.hypot(n[0], n[1]);
      n = [n[0] / l, n[1] / l];
      m = 1 / (n[0] * n1[0] + n[1] * n1[1]);
    } else {
      const d = a ? dir(a, p) : dir(p, b);
      n = [-d[1], d[0]];
    }
    const wl = h * (opts.widthJitter ? 1 + opts.widthJitter * (rnd() - 0.5) : 1);
    const wr = h * (opts.widthJitter ? 1 + opts.widthJitter * (rnd() - 0.5) : 1);
    left.push([p[0] + n[0] * wl * m, p[1] + n[1] * wl * m]);
    right.push([p[0] - n[0] * wr * m, p[1] - n[1] * wr * m]);
  }
  let ring = [...right, ...left.reverse()];
  if (opts.nearPoints) {
    // insert points closer than 0.1 m to an existing vertex (along the edge)
    const out = [];
    for (let i = 0; i < ring.length; i++) {
      out.push(ring[i]);
      if (i % opts.nearPoints === 3) {
        const q = ring[(i + 1) % ring.length];
        const dx = q[0] - ring[i][0], dy = q[1] - ring[i][1];
        const l = Math.hypot(dx, dy);
        const f = (0.03 + 0.05 * rnd()) / l;
        // slightly off the edge so it is not exactly collinear
        out.push([ring[i][0] + dx * f + 0.002 * (rnd() - 0.5), ring[i][1] + dy * f + 0.002 * (rnd() - 0.5)]);
      }
    }
    ring = out;
  }
  if (opts.spikes) {
    // sharp spikes sticking out of the right side (outer boundary), every `spikes` vertices
    const out = [];
    for (let i = 0; i < ring.length; i++) {
      out.push(ring[i]);
      if (i < right.length - 2 && i % opts.spikes === 1) {
        const p = ring[i], q = ring[i + 1];
        const dx = q[0] - p[0], dy = q[1] - p[1];
        const l = Math.hypot(dx, dy);
        const nx = dy / l, ny = -dx / l; // right normal = outward for the right side
        const len = 1.5 + 2 * rnd();
        const base = 0.15 + 0.1 * rnd();
        const mx = p[0] + dx * 0.5, my = p[1] + dy * 0.5;
        out.push([mx - (dx / l) * base, my - (dy / l) * base]);
        out.push([mx + nx * len, my + ny * len]);
        out.push([mx + (dx / l) * base, my + (dy / l) * base]);
      }
    }
    ring = out;
  }
  return ring;
}

// A strip of exactly constant width: both sides offset along the chord normals of a random
// centreline, so the sides are nearly anti-parallel everywhere and collapse at once.
function constStrip(segs, seed, segLen = 2, half = 1.5, curv = 0.03) {
  const rnd = mulberry32(seed);
  const c = [[0, 0]];
  let ang = 0, dang = 0;
  for (let i = 0; i < segs; i++) {
    dang = 0.7 * dang + curv * (rnd() - 0.5);
    ang += dang;
    const [x, y] = c[c.length - 1];
    c.push([x + segLen * Math.cos(ang), y + segLen * Math.sin(ang)]);
  }
  const left = [], right = [];
  for (let i = 0; i < c.length; i++) {
    const a = c[Math.max(i - 1, 0)], b = c[Math.min(i + 1, c.length - 1)];
    const dx = b[0] - a[0], dy = b[1] - a[1], l = Math.hypot(dx, dy);
    const n = [(-dy / l) * half, (dx / l) * half];
    left.push([c[i][0] + n[0], c[i][1] + n[1]]);
    right.push([c[i][0] - n[0], c[i][1] - n[1]]);
  }
  return [...right, ...left.reverse()];
}

function comb(teeth, toothW, gapW, toothH, baseH) {
  // base along x from 0..W, teeth pointing up
  const W = teeth * toothW + (teeth + 1) * gapW;
  const ring = [[0, 0], [W, 0], [W, baseH]];
  let x = W - gapW;
  for (let i = 0; i < teeth; i++) {
    ring.push([x, baseH]);
    ring.push([x, baseH + toothH]);
    ring.push([x - toothW, baseH + toothH]);
    ring.push([x - toothW, baseH]);
    x -= toothW + gapW;
  }
  ring.push([0, baseH]);
  return ring;
}

export function allCases() {
  const cases = [];
  const add = (group, name, rings, note) => cases.push({ group, name, rings, ...(note ? { note } : {}) });

  // --- rectangles and squares ---
  add('rect', 'rect_4x2', [rect(0, 0, 4, 2)]);
  add('rect', 'square_2', [[[0, 0], [2, 0], [2, 2], [0, 2]]], 'four edges meet in one point');
  add('rect', 'square_100', [rect(0, 0, 100, 100)]);
  add('rect', 'rect_10x1', [rect(0, 0, 10, 1)]);
  add('rect', 'rect_3x500', [rect(0, 0, 500, 3)], 'road-like aspect ratio');
  add('rect', 'rect_200x3', [rect(0, 0, 200, 3)]);
  add('rect', 'square_rot30', [rotate(rect(0, 0, 2, 2), Math.PI / 6)]);
  add('rect', 'rect_rot17', [rotate(rect(0, 0, 4, 2), 0.3)]);
  add('rect', 'rect_rot45', [rotate(rect(0, 0, 6, 1), Math.PI / 4)]);
  add('rect', 'rect_cw', [rect(0, 0, 4, 2).reverse()], 'clockwise outer ring; the crate must flip it');
  add('rect', 'rect_closed', [[...rect(0, 0, 4, 2), [0, 0]]], 'closing point present');
  add('rect', 'rect_dup', [[[0, 0], [4, 0], [4, 0], [4, 2], [0, 2], [0, 2]]], 'consecutive duplicate points');
  add('rect', 'rect_collinear', [[[0, 0], [2, 0], [4, 0], [4, 2], [0, 2]]], 'one collinear vertex');
  add('rect', 'rect_collinear_many', [[[0, 0], [1, 0], [2.5, 0], [4, 0], [4, 1], [4, 2], [3, 2], [0, 2], [0, 1]]], 'several collinear vertices');

  // --- triangles and regular polygons ---
  add('regular', 'triangle_4x3', [[[0, 0], [4, 0], [2, 3]]]);
  add('regular', 'triangle_right', [[[0, 0], [4, 0], [0, 3]]]);
  add('regular', 'triangle_thin', [[[0, 0], [10, 0], [5, 0.2]]]);
  for (const n of [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 24, 32, 48, 64]) {
    add('regular', `regular_${n}`, [regular(n, 10, 0.1)]);
  }
  add('regular', 'regular_8_axis', [regular(8, 10, Math.PI / 8)], 'octagon with axis-aligned edges');
  add('regular', 'regular_6_r1000', [regular(6, 1000, 0)]);

  // --- letters ---
  add('letters', 'L', [[[0, 0], [4, 0], [4, 2], [2, 2], [2, 4], [0, 4]]]);
  add('letters', 'L_thin', [[[0, 0], [40, 0], [40, 3], [3, 3], [3, 30], [0, 30]]]);
  add('letters', 'T', [[[2, 0], [4, 0], [4, 4], [6, 4], [6, 6], [0, 6], [0, 4], [2, 4]]], 'bar on top');
  add('letters', 'T_spec', [[[0, 0], [6, 0], [6, 2], [4, 2], [4, 6], [2, 6], [2, 2], [0, 2]]], 'stem below the bar (spec table #4)');
  add('letters', 'U', [[[0, 0], [6, 0], [6, 5], [4, 5], [4, 2], [2, 2], [2, 5], [0, 5]]]);
  add('letters', 'H', [[[0, 0], [2, 0], [2, 3], [5, 3], [5, 0], [7, 0], [7, 8], [5, 8], [5, 4], [2, 4], [2, 8], [0, 8]]]);
  add('letters', 'cross_equal', [[[1, 0], [2, 0], [2, 1], [3, 1], [3, 2], [2, 2], [2, 3], [1, 3], [1, 2], [0, 2], [0, 1], [1, 1]]], 'plus sign, arm width equals arm length');
  add('letters', 'cross_long', [[[4, 0], [6, 0], [6, 4], [10, 4], [10, 6], [6, 6], [6, 10], [4, 10], [4, 6], [0, 6], [0, 4], [4, 4]]]);
  add('letters', 'cross_uneven', [[[4, 0], [5, 0], [5, 4], [10, 4], [10, 7], [5, 7], [5, 10], [4, 10], [4, 7], [0, 7], [0, 4], [4, 4]]]);
  add('letters', 'Z', [[[0, 0], [8, 0], [8, 2], [3, 2], [8, 7], [8, 9], [0, 9], [0, 7], [5, 7], [0, 2]]]);
  add('letters', 'stairs', [[[0, 0], [5, 0], [5, 1], [4, 1], [4, 2], [3, 2], [3, 3], [2, 3], [2, 4], [1, 4], [1, 5], [0, 5]]]);

  // --- combs ---
  add('comb', 'comb_4', [comb(4, 1, 1, 3, 1)], 'teeth as wide as the gaps');
  add('comb', 'comb_8_thin', [comb(8, 0.5, 1.5, 5, 2)]);
  add('comb', 'comb_6_wide', [comb(6, 2, 0.5, 2, 3)]);
  add('comb', 'comb_20', [comb(20, 1, 1, 4, 2)]);
  add('comb', 'comb_double', [[[0, 0], [9, 0], [9, 3], [8, 3], [8, 1], [7, 1], [7, 3], [6, 3], [6, 1], [5, 1], [5, 3], [4, 3], [4, 1], [3, 1], [3, 3], [2, 3], [2, 1], [1, 1], [1, 3], [0, 3]]]);

  // --- stars ---
  add('star', 'star_5', [star(5, 10, 4)]);
  add('star', 'star_7', [star(7, 10, 6, 0.2)]);
  add('star', 'star_12_sharp', [star(12, 10, 2)]);
  add('star', 'star_24', [star(24, 10, 8.5, 0.05)]);
  add('star', 'star_3_spiky', [star(3, 10, 1)]);

  // --- nearly parallel edges ---
  add('parallel', 'trapezoid_taper_1e-3', [[[0, 0], [100, 0], [100, 3.1], [0, 3]]]);
  add('parallel', 'trapezoid_taper_1e-5', [[[0, 0], [500, 0], [500, 3.0005], [0, 3]]]);
  add('parallel', 'wedge_thin', [[[0, 0], [100, -0.5], [100, 0.5]]], 'very sharp convex vertex');
  add('parallel', 'parallelogram_near', [[[0, 0], [100, 0.001], [100, 3.001], [0, 3]]]);
  add('parallel', 'zigzag_strip', [(() => { const r = []; for (let i = 0; i <= 20; i++) r.push([i * 10, (i % 2) * 0.05]); for (let i = 20; i >= 0; i--) r.push([i * 10, 3 + ((i + 1) % 2) * 0.05]); return r; })()]);
  add('parallel', 'bent_strip', [[[0, 0], [100, 0], [200, 1], [200, 4], [100, 3], [0, 3]]]);
  add('parallel', 'almost_collinear', [[[0, 0], [50, 1e-4], [100, 0], [100, 10], [0, 10]]]);

  // --- roads ---
  add('road', 'road_straight_500', [road(50, 1, 10, 3, 0)], 'straight 3 m x 500 m buffered polyline (collinear vertices)');
  add('road', 'road_gentle_50', [road(50, 2, 10, 3, 0.02)]);
  add('road', 'road_curvy_80', [road(80, 3, 5, 4, 0.15)]);
  add('road', 'road_width5_100', [road(100, 4, 6, 5, 0.06, { widthJitter: 0.1 })]);
  add('road', 'road_near_points', [road(40, 5, 10, 3, 0.05, { nearPoints: 7 })], 'points closer than 0.1 m');
  add('road', 'road_spikes', [road(40, 6, 10, 3, 0.04, { spikes: 6 })], 'sharp spikes on one side');
  add('road', 'road_2m_200', [road(200, 7, 3, 2, 0.05, { widthJitter: 0.05 })]);

  add('road', 'const_strip_50', [constStrip(24, 81, 2, 1.5, 0.03)], 'strip of exactly constant width (chord normals)');
  add('road', 'const_strip_100', [constStrip(49, 82, 2, 1.5, 0.1)], 'strip of exactly constant width (chord normals)');
  add('road', 'const_strip_200', [constStrip(99, 83, 2, 2.5, 0.03)], 'strip of exactly constant width (chord normals)');
  add('road', 'const_strip_400', [constStrip(199, 84, 2, 1.0, 0.01)], 'strip of exactly constant width (chord normals)');
  add('bigcoord', 'big_const_strip_300', [translate(constStrip(149, 85, 2, 1.5, 0.05), -30000, 120000)], 'constant-width strip far from the origin');

  // --- holes ---
  add('holes', 'hole_center', [rect(0, 0, 10, 10), [[3, 3], [3, 7], [7, 7], [7, 3]]], 'spec table #6');
  add('holes', 'hole_offset', [rect(0, 0, 10, 10), [[2, 3], [2, 6], [5, 6], [5, 3]]]);
  add('holes', 'hole_ccw_given', [rect(0, 0, 10, 10), [[3, 3], [7, 3], [7, 7], [3, 7]]], 'hole given counter-clockwise; the crate must flip it');
  add('holes', 'hole_near_outer', [rect(0, 0, 10, 10), [[0.05, 2], [0.05, 8], [5, 8], [5, 2]]], 'hole 0.05 from the outer ring');
  add('holes', 'hole_near_outer_2', [rect(0, 0, 20, 6), [[1, 0.2], [1, 5.8], [19, 5.8], [19, 0.2]]], 'thin frame');
  add('holes', 'hole_triangle', [regular(6, 10), regular(3, 4, 0.3).reverse()]);
  add('holes', 'hole_rotated', [rect(-10, -10, 20, 20), rotate(rect(-3, -3, 6, 6), 0.4).reverse()]);
  add('holes', 'holes_2', [rect(0, 0, 20, 10), rect(2, 2, 6, 6).reverse(), rect(12, 2, 6, 6).reverse()]);
  add('holes', 'holes_3_grid', [rect(0, 0, 30, 10), rect(2, 2, 6, 6).reverse(), rect(12, 2, 6, 6).reverse(), rect(22, 2, 6, 6).reverse()]);
  add('holes', 'holes_4_random', [rect(0, 0, 40, 40), ...[[5, 5], [25, 8], [8, 27], [27, 26]].map(([x, y], i) => translate(randomStar(9, 100 + i, 0.4, 1, 4), x + 4, y + 4).reverse())]);
  add('holes', 'hole_in_star', [randomStar(30, 7, 0.7, 1, 20), regular(12, 4, 0.1).reverse()]);
  add('holes', 'hole_in_L', [[[0, 0], [10, 0], [10, 4], [4, 4], [4, 10], [0, 10]], [[1, 1], [1, 3], [3, 3], [3, 1]]]);
  add('holes', 'holes_many_small', [rect(0, 0, 50, 20), ...(() => { const hs = []; for (let i = 0; i < 6; i++) for (let j = 0; j < 2; j++) hs.push(regular(5, 1.5, 0.2 * i + j, 5 + 8 * i, 5 + 10 * j).reverse()); return hs; })()]);

  // --- rings touching at a vertex (weakly simple) ---
  add('touch', 'touch_hole_outer', [[[0, 0], [5, 0], [10, 0], [10, 10], [0, 10]], [[5, 0], [3, 4], [7, 4]]], 'triangular hole touching the outer ring at (5, 0)');
  add('touch', 'touch_hole_corner', [rect(0, 0, 10, 10), [[0, 0], [2, 5], [5, 2]].reverse()], 'hole touching an outer corner');
  add('touch', 'touch_two_holes', [rect(0, 0, 20, 10), [[2, 2], [2, 8], [10, 5]].reverse().reverse(), [[10, 5], [18, 8], [18, 2]]], 'two holes sharing a vertex');
  add('touch', 'touch_three_holes', [rect(0, 0, 20, 20), [[10, 10], [4, 14], [4, 6]], [[10, 10], [16, 6], [16, 14]], [[10, 10], [6, 16], [14, 16]]], 'three holes meeting at one point');

  // --- random simple polygons ---
  for (const [n, seed] of [[8, 11], [10, 12], [12, 13], [15, 14], [20, 15], [25, 16], [30, 17], [40, 18], [50, 19], [60, 20], [80, 21], [100, 22], [150, 23], [200, 24]]) {
    add('random', `random_star_${n}_s${seed}`, [randomStar(n, seed)]);
  }
  for (const [n, seed] of [[8, 31], [10, 32], [12, 33], [15, 34], [20, 35], [25, 36], [30, 37], [40, 38], [50, 39], [60, 40], [80, 41], [100, 42]]) {
    add('random', `random_2opt_${n}_s${seed}`, [random2opt(n, seed)]);
  }
  for (const [cells, seed] of [[6, 51], [10, 52], [15, 53], [20, 54], [30, 55], [40, 56], [60, 57], [80, 58]]) {
    add('random', `polyomino_${cells}_s${seed}`, [randomPolyomino(cells, seed)], 'random rectilinear polygon on an integer grid');
  }

  // --- large coordinates (plane rectangular coordinate system, around x=-30000, y=120000) ---
  const BX = -30000, BY = 120000;
  add('bigcoord', 'big_rect', [translate(rect(0, 0, 40, 25), BX, BY)]);
  add('bigcoord', 'big_road_500', [translate(rotate(rect(0, 0, 500, 3), 0.37), BX + 123.4, BY - 56.7)]);
  add('bigcoord', 'big_road_curvy', [translate(road(60, 8, 8, 3.5, 0.08), BX, BY)]);
  add('bigcoord', 'big_L', [translate([[0, 0], [40, 0], [40, 20], [20, 20], [20, 40], [0, 40]], BX, BY)]);
  add('bigcoord', 'big_hole', [translate(rect(0, 0, 100, 100), BX, BY), translate([[30, 30], [30, 70], [70, 70], [70, 30]], BX, BY)]);
  add('bigcoord', 'big_random_60', [translate(randomStar(60, 61, 0.3, 1, 80), BX, BY)]);
  add('bigcoord', 'big_spikes', [translate(road(30, 9, 10, 3, 0.04, { spikes: 5, nearPoints: 9 }), BX, BY)]);

  // --- large vertex counts (CGAL cannot solve these; checked by properties) ---
  add('large', 'ellipse_500', [noisyCircle(500, 71, 100, 0.02, 1.3)]);
  add('large', 'noisy_1000', [noisyCircle(1000, 72, 100, 0.02)]);
  add('large', 'random_star_1000', [randomStar(1000, 73, 0.5, 1, 100)]);
  add('large', 'noisy_2000', [noisyCircle(2000, 74, 100, 0.02)]);
  add('large', 'road_2000', [road(1000, 75, 2, 3, 0.03, { widthJitter: 0.05 })], 'road strip with 2,002 vertices');
  add('large', 'noisy_5000', [noisyCircle(5000, 76, 100, 0.02)]);
  add('large', 'random_star_5000', [randomStar(5000, 77, 0.5, 1, 100)]);

  // --- real data: road-network faces from OpenStreetMap (see oracle/inputs/README.md)
  for (const f of readdirSync(inputsDir).filter((f) => f.endsWith('.json')).sort()) {
    const rings = JSON.parse(readFileSync(join(inputsDir, f), 'utf8'));
    add('osm', f.replace(/\.json$/, ''), rings, 'OpenStreetMap road-network face (ODbL), plane rectangular zone VII');
  }

  return cases;
}
