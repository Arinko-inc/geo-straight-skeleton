# Real-data inputs

`osm_area_*.json` are road-network faces (rings as `[[[x, y], ...], ...]`, ring 0 outer and
counter-clockwise, the others holes and clockwise, no closing point) built from OpenStreetMap
roads, in the Japan Plane Rectangular Coordinate System zone VII (metres). Several have holes
that touch the outer ring or another hole at a single vertex.

Data © OpenStreetMap contributors, available under the Open Database License (ODbL 1.0),
https://www.openstreetmap.org/copyright. These files and the corpus cases built from them
(`tests/corpus/osm_*.json`) are test fixtures only and are not part of the published crate.
