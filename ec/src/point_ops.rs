//! Generic point abstraction.
//!
//! Points are parameterized by their curve model through the associated type
//! `Curve`. This lets us reuse the same interface across Weierstrass,
//! Montgomery, Edwards, etc., while keeping model-specific formulas inside
//! each implementation.
use fp::field_ops::FieldOps;
use subtle::{ConditionallySelectable};

/// Generic group interface for curve points.
///
/// We intentionally do **not** require the standard operator traits here
/// (`Add`, `Sub`, `Mul`, `Neg`) because point addition and negation usually
/// need access to the curve parameters. The clean abstraction boundary is a
/// method-based API taking `&Self::Curve` explicitly.
pub trait PointOps: Clone  + ConditionallySelectable{
    type BaseField: FieldOps;
    type Curve;

    fn identity(curve: &Self::Curve) -> Self;
    fn is_identity(&self) -> bool;
    fn negate(&self, curve: &Self::Curve) -> Self;

    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self;
}

/// Extension trait for points that support full group addition.
///
/// Not every point representation can add two arbitrary points (e.g.
/// Montgomery x-only points). Protocols that need addition (like ElGamal)
/// should bound on `PointAdd` instead of plain `PointOps`.
pub trait PointAdd: PointOps {
    fn add(&self, other: &Self, curve: &Self::Curve) -> Self;
}

