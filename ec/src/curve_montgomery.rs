//! Elliptic curve definition (Montgomery form).
//!
//! # Equation
//!
//! A Montgomery curve over a field `F` is given by
//!
//! ```text
//! B y² = x(x² + A x + 1)
//! ```
//!
//! where `B ≠ 0`.
//!
//! # Representation choice
//!
//! In this module, the curve parameters are stored as the pair `(A, B)`.
//! The native point representation for arithmetic is **x-only projective**
//! coordinates on the Kummer line, rather than full affine points.
//!
//! # Why x-only arithmetic?
//!
//! On a Montgomery curve, scalar multiplication can be implemented using only
//! x-coordinates via the Montgomery ladder. This works on the Kummer quotient
//! `E / {±1}`, where a point `P` and its inverse `-P` have the same image.
//!
//! The advantage is that the ladder uses a uniform sequence of differential
//! additions and doublings, which is especially convenient for constant-time
//! scalar multiplication.


use fp::field_ops::FieldOps;

use crate::curve_ops::Curve;
use crate::point_montgomery::KummerPoint;

/// A Montgomery curve
///
/// ```text
/// B y² = x(x² + A x + 1)
/// ```
///
/// over a field `F`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MontgomeryCurve<F: FieldOps> {
    pub a: F,
    pub b: F,
}

impl<F: FieldOps> MontgomeryCurve<F> {
    // -------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------

    /// Construct a Montgomery curve from its two coefficients `A` and `B`.
    pub fn new(a: F, b: F) -> Self {
        assert!(Self::montgomery_is_smooth(&a, &b));
        Self { a, b }
    }

    pub fn montgomery_is_smooth(a: &F, b: &F) -> bool {
        let two = F::one().double();
        !bool::from(b.is_zero()) && *a != two && *a != two.negate()
    }

    // -------------------------------------------------------------------
    // Curve predicates
    // -------------------------------------------------------------------


    // -------------------------------------------------------------------
    // Constants used by x-only formulas
    // -------------------------------------------------------------------

    /// Return the Montgomery-ladder constant traditionally denoted `A24`.
    ///
    /// Depending on the exact doubling formula you choose, this is usually one
    /// of the constants derived from `A`, such as `(A + 2)/4` or `(A - 2)/4`.
    ///
    /// The exact convention should match the doubling formula implemented in
    /// [`KummerPoint::double`].
    pub fn a24(&self) -> F {
        todo!()
    }

    /// Return the model coefficients `[A, B]`.
    pub fn a_invariants(&self) -> [F; 2] {
        todo!()
    }
}

impl<F: FieldOps> Curve for MontgomeryCurve<F> {
    type BaseField = F;
    type Point = KummerPoint<F>;

    /// Return `true` if `point` is a valid Kummer/x-line point for this curve.
    ///
    /// Since the Kummer line forgets the sign of `y`, this is not the same as
    /// checking whether a full affine point lies on the curve.
    ///
    /// A simple policy is to accept every nonzero projective pair `(X:Z)` as a
    /// valid x-line point, with `(1:0)` representing the identity image.
    fn is_on_curve(&self, point: &Self::Point) -> bool {
        todo!()
    }

    fn random_point(&self) -> Self::Point {
        todo!()
    }

    fn j_invariant(&self) -> u64 {
        todo!()
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        todo!()
    }
}