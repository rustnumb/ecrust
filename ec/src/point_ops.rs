//! Generic point abstraction.
//!
//! Points are parameterized by their curve model through the associated type
//! `Curve`. This lets us reuse the same interface across Weierstrass,
//! Montgomery, Edwards, etc., while keeping model-specific formulas inside
//! each implementation.
use fp::field_ops::FieldOps;
use subtle::ConditionallySelectable;
use crate::point_weierstrass::AffinePoint;

/// Generic group interface for curve points.
///
/// We intentionally do **not** require the standard operator traits here
/// (`Add`, `Sub`, `Mul`, `Neg`) because point addition and negation usually
/// need access to the curve parameters. The clean abstraction boundary is a
/// method-based API taking `&Self::Curve` explicitly.
pub trait PointOps: Clone {
    type BaseField: FieldOps;
    type Curve;

    fn identity(curve: &Self::Curve) -> Self;
    fn is_identity(&self) -> bool;
    fn negate(&self, curve: &Self::Curve) -> Self;
    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self;
}


impl<F> Default for AffinePoint<F>
where
    F: FieldOps + Copy,
{
    fn default() -> Self {Self::identity()}
}