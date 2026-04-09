//! Kernel subgroup of an isogeny.

use ec::point_ops::PointOps;
use ec::curve_ops::Curve;

pub struct KernelSubgroup<C: Curve> {
    pub points: Vec<C::Point>,
}

impl<C: Curve> KernelSubgroup<C> {
    pub fn trivial(curve: &C) -> Self {
        Self { points: vec![C::Point::identity(curve)] }
    }
}