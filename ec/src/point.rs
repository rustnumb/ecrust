//! Elliptic curve point representation and group law.

/// An affine point on an elliptic curve over Fp.
pub struct AffinePoint {
    // TODO: replace with fp::FieldElement once defined
    pub x: u64,
    pub y: u64,
    pub infinity: bool,
}

impl AffinePoint {
    /// Returns the point at infinity (identity element).
    pub fn identity() -> Self {
        Self { x: 0, y: 0, infinity: true }
    }
}
