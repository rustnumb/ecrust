//! Generic curve abstraction.
//!
//! The main goal of this trait is to decouple:
//! 1. the **curve model** (short Weierstrass, Montgomery, Edwards, ...), and
//! 2. the **point representation** (affine, projective, x-only, ...).
//!
//! Each concrete curve type chooses its base field and its native point type
//! through associated types.

use crate::point_ops::PointOps;
use fp::field_ops::FieldOps;

/// Generic elliptic-curve model.
///
/// A curve model fixes:
/// - the base field,
/// - the point type used with that model,
/// - and the membership test `is_on_curve`.
///
/// This is intentionally minimal so that different models (Montgomery,
/// short/general Weierstrass, Edwards, Hessian, ...)
/// can implement it without being forced into one particular formula set.
pub trait Curve: Sized + Clone + PartialEq + Eq {
    /// Base field of the curve.
    type BaseField: FieldOps;

    /// Native point representation for this curve model.
    type Point: PointOps<Curve = Self, BaseField = Self::BaseField>;

    /// Return `true` if `point` is a valid point on this curve.
    fn is_on_curve(&self, point: &Self::Point) -> bool;

    /// Return a random point that is on the curve.
    fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point;

    /// Return the j_invariant of the curve;
    fn j_invariant(&self) -> Self::BaseField;

    // Return the a-invariants of the curve.
    ///
    /// The interpretation depends on the curve model:
    ///   - **Weierstrass** \to `[a_1, a_2, a_3, a_4, a_6]`  (5 elements)
    ///   - **Montgomery**  \to `[A, B]`                   (2 elements)
    ///   - **Edwards**     \to `[a, d]`                   (2 elements)
    fn a_invariants(&self) -> Vec<Self::BaseField>;

    /// Return the group identity.
    fn identity(&self) -> Self::Point {
        Self::Point::identity(self)
    }
}
