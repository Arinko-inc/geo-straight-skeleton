# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the crate follows
[Semantic Versioning](https://semver.org/).

## 0.1.1 - 2026-10-28

- Add the repository link (the GitHub repository is now public).

## [0.1.0] - Unreleased

### Added

- `Skeleton::build`: straight skeleton of a polygon with holes in `f64`, with one face per input
  edge, skeleton nodes merged by position and time, and `SkeletonError` for invalid input
  (no panics).
- Input handling: either ring orientation, optional closing point, consecutive duplicates,
  collinear points (each keeps its face), rings touching at a shared vertex.
- `Skeleton::arcs`: interior arcs as vertex index pairs.
- `Skeleton::offset`: the inward offset polygon at a given distance, as a `Result` (a contour
  that fails to close is an error, not a ring).
- Feature `geo-types` (geo-types 0.7.8 or later): `Skeleton::from_polygon` and
  `Skeleton::arcs_multi_line_string`.
- `examples/svg.rs`: draw a skeleton and its offsets as SVG.
- Tests against CGAL's `Straight_skeleton_2` on 153 inputs (150 compared, all matching), property
  tests on random polygons, and timing tests.
