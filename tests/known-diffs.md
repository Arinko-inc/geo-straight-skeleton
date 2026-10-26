# Known differences from CGAL and open issues

The oracle test (`tests/oracle.rs`) reads lines of the form `- case: <name>` from this file and
treats those cases as known differences. There are none at the moment.

## Inputs CGAL could not solve (checked by `tests/properties.rs` only)

3 of 148 corpus cases. CGAL's WASM build (npm `straight-skeleton` 2.0.1) fails on them:

- `big_road_curvy`: returns null.
- `noisy_5000`: `std::bad_alloc` (the WASM heap is limited to 2 GB).
- `zigzag_strip`: aborts with `Assertion failed: this->is_initialized()` (boost::optional inside CGAL).

The wrapper also returns null (or runs out of memory) for every input whose rings touch at a
vertex (4 synthetic `touch_*` cases and the OpenStreetMap faces `osm_area_1`, `872`, `1276`, `2278`,
`4216`). Those 9 are compared with CGAL's result on a copy whose touching hole vertices are pulled
into their holes by 2^-20 of the frame (`cgal_perturbed` in the corpus), with a tolerance of 10x
that distance.

## Status (work paused 2026-10-27)

Done:
- Crate renamed to `geo-straight-skeleton`; feature `geo-types` (`Skeleton::from_polygon`,
  `Skeleton::arcs_multi_line_string`); `Skeleton::offset`; `examples/svg.rs`; `tests/perf.rs`.
- Oracle: 136 cases match CGAL exactly (face areas, vertex sets, arc counts, per-face vertex
  sequences), 9 more match the perturbed copy. Property tests pass on all 148 cases and on
  50,000 random cases per generator (except the open issue below).
- Rings touching at a vertex are accepted. 5,000 vertices take about 75 ms (release, M5 Max).

Open issue: strips of exactly constant width built from chord normals
(`constant_width_strips` in `tests/properties.rs`, ignored for now). The whole centreline collapses
within ~1e-11 of the same time, and sliver tips (sector ~1e-9 rad, speed ~1e9) get positions that
are off along the sliver because their creation time carries the rounding of an ill-conditioned
event. Symptom: a spurious node next to the end of the strip (e.g. seed 7360304889168602710,
76 segments, half-width 1.0, curvature 0.03), time slightly above the distance to the boundary.
The perf-test strips (2,000 to 5,000 vertices) no longer fail after the sector-consistency check
and the sliver fallback in `Engine::plan`/`surgery`, but the spurious-node case remains.

Next:
1. Fix the open issue. Anchoring new vertices at their two lines' intersection was tried and did
   not help (reverted). Ideas: treat a sliver tip whose lines both pass through the event point as
   a participant (topological rule, like `find_on_line`); or close slivers instantly but give the
   far node the time and position of its own three-line point.
2. Add constant-width strips to the CGAL corpus (CGAL is slow on them: about 2 s at 400 vertices).
3. Publishing: README (English, results table), CHANGELOG, LICENSE-MIT, LICENSE-APACHE,
   `cargo clippy --all-features -- -D warnings`, `cargo doc`, `cargo package --list`.
