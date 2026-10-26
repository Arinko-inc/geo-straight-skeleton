// Builds the oracle corpus: tests/corpus/<name>.json = input rings + CGAL's straight skeleton.
//
// The oracle is the npm package `straight-skeleton` 2.0.1 (CGAL Straight_skeleton_2 compiled to WASM),
// taken from ~/Algo-main/node_modules (override with ALGO_MAIN=<dir>). Usage:
//
//   node oracle/gen.mjs            # (re)generate every case from oracle/shapes.mjs, run CGAL on each
//   node oracle/gen.mjs rect_4x2   # only the named cases
//
// The WASM wrapper passes coordinates as float32 both ways. To keep CGAL's input identical to ours,
// every input point is snapped so that its local coordinate (p - c) / s is a multiple of 2^-24 with
// |local| < 1 (c = centre, s = power of two >= the extent); such values are exact in float32, and
// p = c + s * local is exact in float64. CGAL then sees exactly our polygon, up to the exact scaling.
//
// CGAL's output vertices come back rounded to float32 (about 6e-8 of the extent). We store them as
// `raw`, and also a refined copy `vertices`: every skeleton vertex lies at distance `time` from the
// supporting line of each face it belongs to, so (x, y, time) is recomputed in float64 by least squares
// from those lines (CGAL's topology, float64 geometry). A refined vertex is kept only if it is within
// 1e-5 of the extent from the raw one; otherwise the raw value is kept and `refined[i]` is false.

import { spawnSync } from 'child_process';
import { writeFileSync, mkdirSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';
import { allCases } from './shapes.mjs';

const here = dirname(fileURLToPath(import.meta.url));
const outDir = join(here, '..', 'tests', 'corpus');
mkdirSync(outDir, { recursive: true });

const area = (r) => {
  let a = 0;
  for (let i = 0; i < r.length; i++) {
    const p = r[i], q = r[(i + 1) % r.length];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return a / 2;
};

function cleanRing(r) {
  const out = [];
  for (const p of r) {
    const last = out[out.length - 1];
    if (last && last[0] === p[0] && last[1] === p[1]) continue;
    out.push(p);
  }
  while (out.length > 1 && out[0][0] === out[out.length - 1][0] && out[0][1] === out[out.length - 1][1]) out.pop();
  return out;
}

function frame(rings) {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const r of rings) for (const [x, y] of r) {
    x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y);
  }
  const ext = Math.max(x1 - x0, y1 - y0);
  const s = Math.pow(2, Math.ceil(Math.log2(ext)));
  const q = s / 1024;
  const cx = Math.round((x0 + x1) / 2 / q) * q, cy = Math.round((y0 + y1) / 2 / q) * q;
  return { cx, cy, s, ext, diag: Math.hypot(x1 - x0, y1 - y0) };
}

const G = 2 ** 24;
const toLocal = (f, [x, y]) => [Math.round(((x - f.cx) / f.s) * G) / G, Math.round(((y - f.cy) / f.s) * G) / G];
const toGlobal = (f, [x, y]) => [f.cx + f.s * x, f.cy + f.s * y];

function runCgal(localRings) {
  const res = spawnSync(process.execPath, [join(here, 'run-cgal.mjs')], {
    input: JSON.stringify({ rings: localRings }),
    maxBuffer: 1 << 30,
    timeout: 120000,
    encoding: 'utf8',
  });
  if (res.error && res.error.code === 'ETIMEDOUT') return { status: 'timeout', message: 'no answer within 120 s' };
  if (res.status !== 0 || !res.stdout) {
    const msg = (res.stderr || '').split('\n').find((l) => /Aborted|Error|error/.test(l)) || `exit ${res.status}`;
    return { status: 'error', message: msg.slice(0, 300) };
  }
  return JSON.parse(res.stdout);
}

// Solve the least-squares system for (x, y, t): n_i . p - t = c_i for each incident face line.
// Householder QR on the m x 3 matrix (the normal equations lose too many digits when the lines are
// nearly concurrent-tangent).
function refineVertex(lines) {
  const m = lines.length;
  const A = lines.map(({ n }) => [n[0], n[1], -1]);
  const b = lines.map(({ c }) => c);
  for (let k = 0; k < 3; k++) {
    let norm = 0;
    for (let i = k; i < m; i++) norm += A[i][k] * A[i][k];
    norm = Math.sqrt(norm);
    if (norm < 1e-300) return null;
    const alpha = A[k][k] > 0 ? -norm : norm;
    const v = new Array(m).fill(0);
    for (let i = k; i < m; i++) v[i] = A[i][k];
    v[k] -= alpha;
    let vv = 0;
    for (let i = k; i < m; i++) vv += v[i] * v[i];
    if (vv < 1e-300) continue;
    for (let j = k; j < 3; j++) {
      let s = 0;
      for (let i = k; i < m; i++) s += v[i] * A[i][j];
      const f = (2 * s) / vv;
      for (let i = k; i < m; i++) A[i][j] -= f * v[i];
    }
    let s = 0;
    for (let i = k; i < m; i++) s += v[i] * b[i];
    const f = (2 * s) / vv;
    for (let i = k; i < m; i++) b[i] -= f * v[i];
  }
  // back substitution on the 3x3 upper triangle; reject (near) rank deficiency
  const R = A;
  const scale = Math.max(Math.abs(R[0][0]), Math.abs(R[1][1]), Math.abs(R[2][2]));
  for (let k = 0; k < 3; k++) if (Math.abs(R[k][k]) < 1e-9 * scale) return null;
  const z = [0, 0, 0];
  for (let k = 2; k >= 0; k--) {
    let s = b[k];
    for (let j = k + 1; j < 3; j++) s -= R[k][j] * z[j];
    z[k] = s / R[k][k];
  }
  return z;
}

function solve(local, f) {
  const res = runCgal(local);
  const cg = { status: res.status, ms: res.ms !== undefined ? Math.round(res.ms * 1000) / 1000 : undefined };
  if (res.message) cg.message = res.message;
  if (res.status === 'ok') {
    const flat = local.flat();
    const nIn = flat.length;
    // the line of each face: CGAL lists a face as [end, ...skeleton..., start] of its input edge
    const faceLine = res.polygons.map((poly) => {
      const a = poly[poly.length - 1], b = poly[0];
      if (a < nIn && b < nIn) {
        const p = flat[a], q = flat[b];
        const dx = q[0] - p[0], dy = q[1] - p[1];
        const l = Math.hypot(dx, dy);
        const n = [-dy / l, dx / l];
        return { n, c: n[0] * p[0] + n[1] * p[1] };
      }
      return null;
    });
    const incident = res.vertices.map(() => []);
    res.polygons.forEach((poly, k) => { for (const v of poly) incident[v].push(k); });
    const refined = [];
    const verts = res.vertices.map(([x, y, t], i) => {
      if (i < nIn) { refined.push(true); const g = toGlobal(f, flat[i]); return [g[0], g[1], 0]; }
      const lines = [...new Set(incident[i])].map((k) => faceLine[k]).filter(Boolean);
      const sol = lines.length >= 3 ? refineVertex(lines) : null;
      const tol = 1e-5;
      if (sol && Math.hypot(sol[0] - x, sol[1] - y) <= tol && Math.abs(sol[2] - t) <= tol) {
        refined.push(true);
        const g = toGlobal(f, [sol[0], sol[1]]);
        return [g[0], g[1], sol[2] * f.s];
      }
      refined.push(false);
      const g = toGlobal(f, [x, y]);
      return [g[0], g[1], t * f.s];
    });
    cg.vertices = verts;
    cg.raw = res.vertices.map(([x, y, t]) => { const g = toGlobal(f, [x, y]); return [g[0], g[1], t * f.s]; });
    cg.refined = refined;
    cg.polygons = res.polygons;
    cg.unrefined = refined.filter((r) => !r).length;
  }
  return cg;
}

// Pull every hole vertex that coincides with a vertex of another ring into its hole (along the
// bisector of the hole's angle there) by 2^-20 local units. Returns null when nothing touches.
function perturbTouches(local) {
  const key = (p) => `${p[0]},${p[1]}`;
  const count = new Map();
  local.forEach((r, ri) => r.forEach((p) => {
    const k = key(p);
    if (!count.has(k)) count.set(k, new Set());
    count.get(k).add(ri);
  }));
  const delta = 2 ** -20;
  let moved = 0;
  const rings = local.map((r, ri) => {
    if (ri === 0) return r.slice();
    return r.map((p, i) => {
      if (count.get(key(p)).size < 2) return p;
      const a = r[(i + r.length - 1) % r.length], b = r[(i + 1) % r.length];
      const u1 = [a[0] - p[0], a[1] - p[1]], u2 = [b[0] - p[0], b[1] - p[1]];
      const l1 = Math.hypot(...u1), l2 = Math.hypot(...u2);
      let bis = [u1[0] / l1 + u2[0] / l2, u1[1] / l1 + u2[1] / l2];
      const turn = (p[0] - a[0]) * (b[1] - p[1]) - (p[1] - a[1]) * (b[0] - p[0]);
      let lb = Math.hypot(...bis);
      if (lb < 1e-9) { bis = [b[1] - a[1], -(b[0] - a[0])]; lb = Math.hypot(...bis); }
      // holes are clockwise: a right turn (turn < 0) means the hole's angle there is convex
      const sgn = turn < 0 || lb < 1e-9 ? 1 : -1;
      moved++;
      return [Math.round((p[0] + (sgn * delta * bis[0]) / lb) * G) / G, Math.round((p[1] + (sgn * delta * bis[1]) / lb) * G) / G];
    });
  });
  return moved ? { rings, delta, moved } : null;
}

function processCase(c) {
  const t0 = Date.now();
  const f = frame(c.rings.map(cleanRing));
  // snap the input as given (closing points, duplicates and orientation are kept as the case wrote them)
  const rings = c.rings.map((r) => r.map((p) => toGlobal(f, toLocal(f, p))));
  // CGAL wants: no closing point, no duplicates, outer CCW, holes CW
  const local = rings.map((r) => cleanRing(r).map((p) => toLocal(f, p))).map((r, i) => {
    const a = area(r);
    return (i === 0 ? a < 0 : a > 0) ? r.slice().reverse() : r;
  });
  const out = { name: c.name, group: c.group, ...(c.note ? { note: c.note } : {}), rings };
  if (local.some((r) => r.length < 3)) {
    out.cgal = { status: 'skipped', message: 'ring with fewer than 3 points' };
    return out;
  }
  const cg = solve(local, f);
  out.cgal = cg;
  out.frame = { cx: f.cx, cy: f.cy, s: f.s };
  // Rings touching at a vertex: CGAL's wrapper rejects them, so also solve a copy in which every
  // touching hole vertex is pulled into its hole by `delta` (the tests compare with a tolerance
  // scaled to it).
  const pert = perturbTouches(local);
  if (pert) {
    const pc = solve(pert.rings, f);
    pc.delta = pert.delta * f.s;
    pc.moved = pert.moved;
    out.cgal_perturbed = pc;
    cg.message = (cg.message ? cg.message + '; ' : '') + `perturbed copy: ${pc.status}`;
  }
  console.log(`${c.name.padEnd(28)} ${String(rings.reduce((a, r) => a + r.length, 0)).padStart(5)} pts  cgal=${cg.status}${out.cgal_perturbed ? `/${out.cgal_perturbed.status}` : ''}${cg.unrefined ? ` unrefined=${cg.unrefined}` : ''}${cg.message ? ` (${cg.message.slice(0, 80)})` : ''}  ${Date.now() - t0} ms`);
  return out;
}

const only = new Set(process.argv.slice(2));
for (const c of allCases()) {
  if (only.size && !only.has(c.name)) continue;
  const out = processCase(c);
  writeFileSync(join(outDir, `${c.name}.json`), JSON.stringify(out));
}
