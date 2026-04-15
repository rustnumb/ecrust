//! Isogeny map between two elliptic curves.

use ec::curve_ops::Curve;

/// The wrapper structure for an isogeny of elliptic curves
pub struct Isogeny<C: Curve> {
    /// The domain of the isogeny, an elliptic curve of some sort
    pub domain: C,
    /// The codomain of the isogeny, an elliptic curve of some sort
    pub codomain: C,
    /// The degree of the isogeny
    pub degree: u64,
}

impl<C: Curve> Isogeny<C> {
    /// Evaluate the isogeny at a point
    pub fn evaluate(&self, _p: &C::Point) -> Option<C::Point> {
        todo!("Vélu's formulas not yet implemented")
    }
}
