//! Kernel subgroup of an isogeny.

use ec::curve_ops::Curve;
use ec::point_ops::PointOps;

/// The structure for a kernel of an isogeny
pub struct KernelSubgroup<C: Curve> {
    /// All the points in the kernel of the isogney
    pub points: Vec<C::Point>,
}

impl<C: Curve> KernelSubgroup<C> {
    /// Check whether the isogney kernel is trivial
    pub fn trivial(curve: &C) -> Self {
        Self {
            points: vec![C::Point::identity(curve)],
        }
    }
}
