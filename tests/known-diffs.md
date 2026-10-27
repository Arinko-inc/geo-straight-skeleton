# Known differences from CGAL

The oracle test (`tests/oracle.rs`) reads lines of the form `- case: <name>` from this file and
treats those cases as known differences. There are none: all 150 corpus cases that CGAL can solve
(141 directly, 9 through a perturbed copy) match.

## Inputs CGAL could not solve (checked by `tests/properties.rs` only)

3 of 153 corpus cases. CGAL's WebAssembly build (npm `straight-skeleton` 2.0.1) fails on them:

- `big_road_curvy`: returns null.
- `noisy_5000`: `std::bad_alloc` (the WebAssembly heap is limited to 2 GB).
- `zigzag_strip`: aborts with `Assertion failed: this->is_initialized()` (a boost::optional inside
  CGAL).

## Inputs compared with a perturbed copy

The wrapper returns null (or runs out of memory) for every input whose rings touch at a vertex:
the 4 synthetic `touch_*` cases and the OpenStreetMap faces `osm_area_1`, `osm_area_872`,
`osm_area_1276`, `osm_area_2278` and `osm_area_4216`. `oracle/gen.mjs` also runs CGAL on a copy
whose touching hole vertices are pulled into their holes by 2^-20 of the frame (about 1e-6 of the
extent; stored as `cgal_perturbed`), and the oracle test compares with that, allowing ten times
the displacement.

## Limitation (not a difference from CGAL)

Floating point throughout, no exact predicates. For strips of exactly constant width (both sides
offset along the same chord normals, so the whole centre line collapses within ~1e-11 of one time)
skeleton nodes carry rounding noise of up to about 1e-8 of the extent; this can show as an extra
arc of that length, and `tests/common` checks face monotonicity at 1e-7 of the extent for that
reason. The five constant-width strips in the corpus (50 to 400 vertices) match CGAL.
