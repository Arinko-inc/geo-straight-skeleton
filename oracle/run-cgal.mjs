// Child process: reads {"rings": [[[x,y],...],...]} (already oriented, f32-exact, no closing point)
// from stdin, runs CGAL (npm straight-skeleton 2.0.1) and prints the raw result as JSON.
// Runs in its own process because the WASM build aborts on bad input and can exhaust its 2 GB heap.
import { createRequire } from 'module';
import { readFileSync } from 'fs';

globalThis.self = globalThis;
globalThis.window = globalThis;

const algo = process.env.ALGO_MAIN || `${process.env.HOME}/Algo-main/`;
const require = createRequire(algo.endsWith('/') ? algo : algo + '/');
const { SkeletonBuilder } = require('straight-skeleton');

const input = JSON.parse(readFileSync(0, 'utf8'));
await SkeletonBuilder.init();
const rings = input.rings.map((r) => [...r, r[0]]);
const t0 = performance.now();
let out;
try {
  const res = SkeletonBuilder.buildFromPolygon(rings);
  const ms = performance.now() - t0;
  if (res === null) out = { status: 'null', ms };
  else out = { status: 'ok', ms, vertices: res.vertices, polygons: res.polygons };
} catch (e) {
  out = { status: 'error', ms: performance.now() - t0, message: String(e && e.message ? e.message : e).slice(0, 300) };
}
process.stdout.write(JSON.stringify(out));
