//! Kernel subgroup of an isogeny.

use ec::point::AffinePoint;

/// A finite subgroup of E(Fp) that will be used as the kernel of an isogeny.
pub struct KernelSubgroup {
    pub points: Vec<AffinePoint>,
}

impl KernelSubgroup {
    /// Constructs the trivial kernel (identity only).
    pub fn trivial() -> Self {
        Self { points: vec![AffinePoint::identity()] }
    }
}
