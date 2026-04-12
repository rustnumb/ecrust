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
        assert!(Self::is_smooth(&a, &b));
        Self { a, b }
    }

    pub fn is_smooth(a: &F, b: &F) -> bool {
        if F::characteristic()[0] != 2 {
            let two = F::one().double();
            !bool::from(b.is_zero()) && *a != two && *a != two.negate()
        } else {
            todo!()
        }
    }


    // -------------------------------------------------------------------
    // Invariants and constants
    // -------------------------------------------------------------------

    /// Return the model coefficients `[A, B]`.
    pub fn a_invariants(&self) -> [F; 2] {
        [self.a, self.b]
    }

    /// Return the Montgomery-ladder constant `A24 = (a + 2)/4`.
    pub fn a24(&self) -> F {
        assert_ne!(F::characteristic()[0], 2);
        let two = <F as FieldOps>::double(&F::one());
        let four = <F as FieldOps>::double(&two);
        let four_inv = <F as FieldOps>::invert(&four).unwrap();

        (self.a + two) * four_inv
    }
}


// -------------------------------------------------------------------
// Curve predicates
// -------------------------------------------------------------------

impl<F: FieldOps> Curve for MontgomeryCurve<F> {
    type BaseField = F;
    type Point = KummerPoint<F>;

    /// Return `true` if `point` is a valid Kummer/x-line point for this curve.
    ///
    /// It checks if `(x³ + A x² + x)/B` is a square in F or not.
    fn is_on_curve(&self, point: &Self::Point) -> bool {
        if point.is_identity() {
            return true;
        }
        else {
            if F::characteristic()[0] != 2 {
                // Normalise projective (X : Z) to affine x = X / Z.
                let z_inv = <F as FieldOps>::invert(&point.z);
                if bool::from(z_inv.is_none()) {
                    // Z = 0 means identity, already handled above.
                    return true;
                }
                let x = point.x * z_inv.unwrap();

                let xsq = <F as FieldOps>::square(&x);
                let xcubed = x * xsq;
                let axsq = self.a * xsq;
                let binv = <F as FieldOps>::invert(&self.b).unwrap();
                let sum = xcubed + axsq + x;
                <F as FieldOps>::legendre(&(sum * binv)) > -1i8
            }
            else {
                todo!()
            }
        }
    }

    fn random_point(&self) -> Self::Point {
        todo!()
    }

    fn j_invariant(&self) -> F {
        if F::characteristic()[0] != 2 {
            let two = <F as FieldOps>::double(&F::one());
            let three = two + F::one();
            let four = <F as FieldOps>::double(&two);
            let sixteen = <F as FieldOps>::square(&four);
            let twofivesix = <F as FieldOps>::square(&sixteen);

            let asq = <F as FieldOps>::square(&self.a);
            let asq_min_three = asq - three;
            let asq_min_three_cubed = asq_min_three * <F as FieldOps>::square(&asq_min_three);


            let asq_min_four_inv = <F as FieldOps>::invert(&(asq - four)).into_option().expect("a should be different from 2, -2.");

            twofivesix * asq_min_three_cubed * asq_min_four_inv
        }
        else {
            todo!()
        }
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        MontgomeryCurve::a_invariants(self).to_vec()
    }
}