//! Generic point abstraction.
//!
//! Points are parameterized by their curve model through the associated type
//! `Curve`. This lets us reuse the same interface across Weierstrass,
//! Montgomery, Edwards, etc., while keeping model-specific formulas inside
//! each implementation.
use fp::field_ops::FieldOps;
use subtle::ConditionallySelectable;

/// Generic group interface for curve points.
///
/// We intentionally do **not** require the standard operator traits here
/// (`Add`, `Sub`, `Mul`, `Neg`) because point addition and negation usually
/// need access to the curve parameters. The clean abstraction boundary is a
/// method-based API taking `&Self::Curve` explicitly.
pub trait PointOps: Clone + ConditionallySelectable {
    /// The base field $\mathbb{F}_{p^M}$
    type BaseField: FieldOps;

    /// The elliptic curve we're working on
    type Curve;

    /// Returns the identity
    fn identity(curve: &Self::Curve) -> Self;

    /// Returns true if and only if `self` is the identity
    fn is_identity(&self) -> bool;

    /// Negate a point
    fn negate(&self, curve: &Self::Curve) -> Self;

    /// Scalar multiplication  `[k]P`  (variable-time double-and-add).
    ///
    /// Provided as a default so every `PointOps` implementor gets it
    /// automatically.
    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self;
}

/// Extension trait for points that support full group addition.
///
/// Not every point representation can add two arbitrary points (e.g.
/// Montgomery x-only points). Protocols that need addition (like ElGamal)
/// should bound on `PointAdd` instead of plain `PointOps`.
pub trait PointAdd: PointOps {
    /// Add a pair of points
    fn add(&self, other: &Self, curve: &Self::Curve) -> Self;
}
