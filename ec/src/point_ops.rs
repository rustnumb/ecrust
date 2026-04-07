//! Generic point abstraction.
//!
//! Points are parameterized by their curve model through the associated type
//! `Curve`. This lets us reuse the same interface across Weierstrass,
//! Montgomery, Edwards, etc., while keeping model-specific formulas inside
//! each implementation.
use subtle::{ConditionallySelectable};


/// Generic group interface for curve points.
///
/// We intentionally do **not** require the standard operator traits here
/// (`Add`, `Sub`, `Mul`, `Neg`) because point addition and negation usually
/// need access to the curve parameters. The clean abstraction boundary is a
/// method-based API taking `&Self::Curve` explicitly.
pub trait PointOps: Clone {
    type BaseField;
    type Curve;

    fn identity(curve: &Self::Curve) -> Self;
    fn is_identity(&self) -> bool;
    fn negate(&self, curve: &Self::Curve) -> Self;
    fn add(&self, rhs: &Self, curve: &Self::Curve) -> Self;
    fn double(&self, curve: &Self::Curve) -> Self;

    /// Scalar multiplication  `[k]P`  (variable-time double-and-add).
    ///
    /// Provided as a default so every `PointOps` implementor gets it
    /// automatically.  For constant-time scalar multiplication, see
    /// [`CtPointOps::scalar_mul_ct`].
    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        let mut result = Self::identity(curve);
        let mut found_one = false;

        for &limb in k.iter().rev() {
            for bit in (0..64).rev() {
                if found_one {
                    result = result.double(curve);
                }
                if (limb >> bit) & 1 == 1 {
                    found_one = true;
                    result = result.add(self, curve);
                }
            }
        }

        result
    }
}



pub trait CtPointOps: PointOps + ConditionallySelectable {
    fn scalar_mul_ct(&self, k: &[u64], curve: &Self::Curve) -> Self;
}
