//! Conversions from and to [`geo_types`] (feature `geo-types`).

use crate::{Skeleton, SkeletonError};
use geo_types::{Coord, LineString, MultiLineString, Polygon};

fn ring(ls: &LineString<f64>) -> Vec<[f64; 2]> {
    ls.0.iter().map(|c| [c.x, c.y]).collect()
}

impl Skeleton {
    /// Straight skeleton of a `geo_types` polygon: the exterior is the outer ring, the
    /// interiors are holes. Orientation and the closing coordinate are handled as in
    /// [`Skeleton::build`]; input point `k` is the `k`-th coordinate of the exterior (without the
    /// closing one), followed by the interiors in order.
    ///
    /// ```
    /// use geo_types::polygon;
    /// use geo_straight_skeleton::Skeleton;
    ///
    /// let p = polygon![(x: 0.0, y: 0.0), (x: 4.0, y: 0.0), (x: 4.0, y: 2.0), (x: 0.0, y: 2.0)];
    /// let sk = Skeleton::from_polygon(&p).unwrap();
    /// assert_eq!(sk.arcs_multi_line_string().0.len(), 5);
    /// ```
    pub fn from_polygon(polygon: &Polygon<f64>) -> Result<Skeleton, SkeletonError> {
        let mut rings = Vec::with_capacity(1 + polygon.interiors().len());
        rings.push(ring(polygon.exterior()));
        for h in polygon.interiors() {
            rings.push(ring(h));
        }
        Skeleton::build(&rings)
    }

    /// The interior arcs of the skeleton ([`Skeleton::arcs`]) as two-point line strings.
    pub fn arcs_multi_line_string(&self) -> MultiLineString<f64> {
        let c = |i: usize| Coord { x: self.vertices[i][0], y: self.vertices[i][1] };
        MultiLineString::new(self.arcs().into_iter().map(|(a, b)| LineString::new(vec![c(a), c(b)])).collect())
    }
}
