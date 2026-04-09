//! Isogeny map between two elliptic curves.

use ec::curve_ops::Curve;

pub struct Isogeny<C: Curve> {
    pub domain:    C,
    pub codomain:  C,
    pub degree:    u64,
}

impl<C: Curve> Isogeny<C> {
    pub fn evaluate(&self, _p: &C::Point) -> Option<C::Point> {
        todo!("Vélu's formulas not yet implemented")
    }
}
