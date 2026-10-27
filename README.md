# geo-straight-skeleton

The straight skeleton of polygons with holes, in `f64`, in pure Rust with no dependencies.

The straight skeleton is what the wavefront of a polygon traces when every edge moves inwards
at unit speed: each input edge sweeps one face, and the faces meet along the skeleton's arcs.
It is the "roof" of a building outline (every face slopes up from its edge at 45°), and the
basis of offset polygons and of medial-axis-like centre lines of roads and channels.

```text
 (0,2) +-------------------------+ (4,2)        faces: one per input edge
       | \          top        / |              time:  0 on the boundary,
       |   \                 /   |                     1 on the ridge
       |left +-------------+ right
       |   /  (1,1)   (3,1)  \   |              arcs:  4 from the corners,
       | /        bottom       \ |                     1 ridge (1,1)-(3,1)
 (0,0) +-------------------------+ (4,0)
```

- Robust against simultaneous events: squares and regular polygons collapse into one node,
  rectilinear shapes, combs, crosses, thin road strips with nearly parallel sides, and holes that
  touch the outer ring at a vertex.
- Checked against CGAL's `Straight_skeleton_2` on 153 inputs (see [Verification](#verification)).
- About 75 ms for 5,000 vertices; memory linear in the input size.
- `#![forbid(unsafe_code)]`, no runtime dependencies; optional [`geo-types`] conversions.

## Usage

```rust
use geo_straight_skeleton::Skeleton;

// ring 0 is the outer boundary, the other rings are holes (any orientation, closing point optional)
let outer = vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]];
let hole = vec![[3.0, 3.0], [3.0, 7.0], [7.0, 7.0], [7.0, 3.0]];
let sk = Skeleton::build(&[outer, hole]).unwrap();

// vertices: the 8 input points (time 0), then the skeleton nodes ([x, y, time])
assert_eq!(sk.vertices.len(), 12);
assert!((sk.vertices[8][2] - 1.5).abs() < 1e-12);
// one face per input edge, counter-clockwise, starting at the edge's end point
assert_eq!(sk.faces.len(), 8);
assert_eq!(sk.faces[0].edge, (0, 1));
// interior arcs as vertex index pairs
assert_eq!(sk.arcs().len(), 12);
// the polygon shrunk by 1: an outer ring and a hole ring
assert_eq!(sk.offset(1.0).unwrap().len(), 2);
```

With the `geo-types` feature:

```rust
# #[cfg(feature = "geo-types")]
# {
use geo_straight_skeleton::Skeleton;
use geo_types::polygon;

let p = polygon![(x: 0.0, y: 0.0), (x: 4.0, y: 0.0), (x: 4.0, y: 2.0), (x: 0.0, y: 2.0)];
let sk = Skeleton::from_polygon(&p).unwrap();
let lines: geo_types::MultiLineString<f64> = sk.arcs_multi_line_string();
assert_eq!(lines.0.len(), 5);
# }
```

`examples/svg.rs` reads a polygon from JSON and draws its skeleton (and offsets) as SVG:

```text
cargo run --example svg -- polygon.json --offsets 0.5,1 > skeleton.svg
```

## Input and output

**Input.** `rings[0]` is the outer boundary, the other rings are holes. Orientation is fixed
internally (outer counter-clockwise, holes clockwise), a closing point equal to the first point is
dropped, and consecutive duplicate points are ignored (closer than 1e-12 of the extent).
Collinear points are kept: each gets its own face, separated from its neighbour by an arc
perpendicular to the edge, as in CGAL. Rings may touch each other at a shared vertex (weakly simple
input, as produced by road-network faces), as long as no edges cross and the wedges at the shared
point do not overlap. Errors (`SkeletonError`, never a panic): no rings, fewer than three distinct
points, NaN or infinite coordinates, zero area, edges that cross or touch away from a shared
vertex, holes outside the outer ring or inside another hole.

**Output.** `Skeleton::vertices` holds `[x, y, time]`: first the input points in input order
(rings one after another, closing points removed, `time` = 0), then the skeleton nodes. `time` is
the offset distance at which the wavefront reaches the node (for a convex polygon, the distance to
the boundary). Nodes at the same position and time are merged into one (CGAL can output several
coincident nodes, e.g. six at the centre of a regular octagon). `Skeleton::faces[k]` is the face of
the input edge that starts at point `k`: `edge` gives its two endpoints with the interior on the
left, and `vertices` lists the face counter-clockwise from `edge.1` to `edge.0`, CGAL's order.
A zero-length input edge (a duplicate point) gets the two-vertex face `[edge.1, edge.0]`.
`Skeleton::arcs` returns every face edge with at least one endpoint at `time > 0`, once.
Where input points coincide (rings touching, duplicates), faces use their own edge's indices,
so an arc can appear once from each of the coincident points.

`Skeleton::offset(d)` returns the polygon shrunk by `d`: the cross-section of the roof at height
`d`, as rings without closing point (counter-clockwise outer rings, clockwise holes). Input points
that were treated as one point (duplicates, rings touching at a vertex) stay one point there.

A sanity check guards every result: the faces must add up to the polygon's area, otherwise `build`
returns `SkeletonError::Numerical` rather than a wrong skeleton. Very sharp input corners (down to
the input tolerance) are simulated like any other vertex, not taken for zero-width overlaps.

## Design

**Method.** An event-driven simulation of the wavefront. The wavefront is a set of cycles of
moving vertices; each wavefront edge ("piece") lies on the offset line of one input edge, and each
vertex moves along the bisector of the interior sector between its two pieces. Two kinds of events
are kept in an indexed priority queue with one entry per vertex and per piece (so memory stays
linear): a piece collapsing (its endpoints meet) and a reflex vertex reaching a piece. Event points
and times are computed from the three supporting lines involved, not from accumulated vertex
positions, as CGAL does; this keeps them accurate even for vertices that move very fast.

Every event is resolved by the same local operation at the event point. It collects every live
vertex and every piece that passes through the point at that time, groups the vertices into
chains (the passes of the wavefront through the point), sorts the incoming and outgoing rays of all
passes by angle and joins each incoming ray to the outgoing ray just before it (clockwise). That is
the only reconnection that keeps the wavefront simple. Simultaneous events at one point (four sides
of a square, the centre of a regular polygon, two reflex vertices meeting head-on, a split that
coincides with an edge collapse) are therefore resolved together, whatever order the queue reports
them in. A new vertex whose sector has zero width, two exactly anti-parallel pieces lying on top of
each other as in a rectangle once its short sides are gone, is "zipped": it slides along the overlap
at the same instant to the end of the shorter piece, where the operation runs again. That is how
ridges between parallel edges arise.

**Why this method.** The deciding requirement was robustness on the inputs of the first user, a
GIS tool that splits road and channel areas: thin strips several metres wide and hundreds of metres
long, nearly parallel sides, points closer than 0.1 m, sharp spikes, and many simultaneous events
in rectilinear and symmetric shapes. Felkel and Obdržálek's algorithm and its variants handle events
one at a time and are known to fail when events coincide. A kinetic triangulation (Aichholzer and
Held) needs a constrained triangulation and flip events, more code with its own degeneracies. CGAL
gets its robustness from exact predicates, which this crate does not have (no dependencies). The
local reconnection above does not depend on the order of coincident events, and only needs angles
of input edge directions, which are exact.

**Numerics.** Coordinates are translated to the centre of the bounding box and scaled so that its
larger side is 1, then mapped back. Things closer than 1e-12 (positions and times, normalised) are
the same; a reconnection that comes out inconsistent (a chain turning by more than half a turn, a
new sector that does not fit inside both of the sectors it is cut from, unmatched rays) means that
rounding put a participant slightly too far away, and the collection radius is widened tenfold up
to 1e-7. A few further rules close slivers between nearly anti-parallel pieces topologically rather
than numerically (their tips can move a billion times faster than the wavefront): a collapsing chain
whose two pieces end at the same vertex takes that vertex with it, and a sliver whose two pieces
already overlap by a rounding error is closed as if it had zero width.

**Complexity.** Each split candidate check is a few flops on contiguous arrays; a changed piece is
checked against every reflex vertex and a new reflex vertex against every piece, so the worst case
is quadratic, in practice about 75 ms for 5,000 vertices (Apple M5 Max, release build).

## Verification

`tests/oracle.rs` compares the result with CGAL's `Straight_skeleton_2`, through the npm package
`straight-skeleton` 2.0.1 (CGAL compiled to WebAssembly), on a corpus generated by
`oracle/gen.mjs`. The wrapper passes coordinates as float32, so inputs are snapped to values that
are exact in float32 after centring and scaling, and CGAL's float32 output nodes are refined in
float64 from the lines of the faces they belong to. For every case, per face (matched by its input
edge):
the area agrees within a relative 1e-6 and the sequence of (merged) vertices is the same; every
node has a counterpart within 1e-6 of the bounding-box diagonal in x, y and time; the arcs and the
total area agree.

| Inputs                                                   | Cases | Compared with CGAL       | Result |
|----------------------------------------------------------|------:|--------------------------|--------|
| rectangles, squares, collinear and duplicate points      |    14 | 14                       | match  |
| triangles, regular polygons with 3 to 64 sides           |    25 | 25                       | match  |
| L, T, U, H, Z, crosses, stairs                           |    11 | 11                       | match  |
| combs, stars                                             |    10 | 10                       | match  |
| nearly parallel edges, thin wedges                       |     7 | 6                        | match  |
| road-like strips (3 m x 500 m, constant width, spikes)   |    11 | 11                       | match  |
| holes (centred, near the outer ring, several)            |    13 | 13                       | match  |
| rings touching at a vertex                               |     4 | 4 (perturbed copy)       | match  |
| random simple polygons (star-shaped, 2-opt, rectilinear) |    34 | 34                       | match  |
| coordinates around (-30000, 120000)                      |     8 | 7                        | match  |
| 500 to 5,000 vertices                                    |     7 | 6                        | match  |
| real road-network faces (OpenStreetMap), up to 2,958 pts |     9 | 4, 5 (perturbed copy)    | match  |
| **total**                                                | **153** | **141, 9 perturbed**   | **150 of 150** |

CGAL cannot solve 3 cases (it returns null, runs out of WebAssembly memory, or fails an internal
assertion); its wrapper also rejects every input whose rings touch at a vertex. Those touching
inputs are compared with CGAL's result on a copy whose touching hole vertices are pulled into
their holes by about 1e-6 of the extent, with a tolerance of ten times that distance. See
`tests/known-diffs.md`.

`tests/properties.rs` checks every corpus case and tens of thousands of random polygons (star
polygons, rectilinear polygons with holes, road strips of varying and of exactly constant width,
holes nearly or exactly touching the outer ring, large coordinates): one face per edge; faces simple,
counter-clockwise, monotone along their edge, covering the polygon without overlap; every node at
distance `time` from the line of each face around it and no farther from the boundary than
`time`; offsets with the area of the faces above them. `tests/spec.rs` holds fixed expectations.
The corpus and the oracle are in the repository, not in the published package.

## Limitations

- Floating point throughout, no exact predicates. For pathological inputs (strips of exactly
  constant width, where the whole centre line collapses within 1e-11 of one time) nodes carry
  rounding noise of up to about 1e-8 of the extent, which can show as an extra arc that short.
- Only unweighted skeletons (all edges move at the same speed). No 3D roof output beyond `time`.
- Rings may touch at shared vertices, not where a vertex lies on another ring's edge.
- Quadratic worst case for the event search (see Complexity).

## Test data

Part of the test corpus (`oracle/inputs/osm_*.json`, `tests/corpus/osm_*.json`) is derived from
OpenStreetMap: data © OpenStreetMap contributors, available under the Open Database License
(<https://www.openstreetmap.org/copyright>). These files are test fixtures in the repository and are
not part of the published crate.

## License

Licensed under either of the Apache License, Version 2.0 (`LICENSE-APACHE`) or the MIT license
(`LICENSE-MIT`), at your option.

[`geo-types`]: https://crates.io/crates/geo-types
