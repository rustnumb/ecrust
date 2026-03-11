//! Isogeny map between two elliptic curves.

use ec::point::AffinePoint;
use ec::curve::WeierstrassCurve;

/// An isogeny φ : E → E′ of a given degree.
pub struct Isogeny {
    pub domain:    WeierstrassCurve,
    pub codomain:  WeierstrassCurve,
    pub degree:    u64,
}

impl Isogeny {
    /// Evaluate the isogeny on a point of the domain curve.
    ///
    /// Returns `None` if `p` is in the kernel (maps to the point at infinity).
    pub fn evaluate(&self, _p: &AffinePoint) -> Option<AffinePoint> {
        // TODO: implement Vélu's formulas
        todo!("Vélu's formulas not yet implemented")
    }
}
