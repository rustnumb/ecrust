//! Elliptic curve definition (general Weierstrass form).
//!
//! # Equation
//!
//! The **general Weierstrass** equation over a field `F` is
//!
//! ```text
//! y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆
//! ```
//!
//! This form is valid in *any* characteristic, including characteristic 2.
//!
//! # Short Weierstrass specialisation
//!
//! When `char(F) ≠ 2, 3` the curve can be brought to the simpler
//!
//! ```text
//! y² = x³ + ax + b          (a₁ = a₂ = a₃ = 0,  a₄ = a,  a₆ = b)
//! ```
//!
//! via the convenience constructor [`WeierstrassCurve::new_short`].

use fp::field_ops::FieldOps;

use crate::curve_ops::Curve;
use crate::point_weierstrass::AffinePoint;

/// An elliptic curve in general Weierstrass form over a field `F`.
///
/// ```text
/// y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆
/// ```
///
/// All five coefficients are stored explicitly; the short Weierstrass case
/// simply has `a1 = a2 = a3 = 0`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeierstrassCurve<F: FieldOps> {
    pub a1: F,
    pub a2: F,
    pub a3: F,
    pub a4: F,
    pub a6: F,
}

impl<F: FieldOps> WeierstrassCurve<F> {
    // -------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------

    /// Construct a curve from all five general Weierstrass coefficients.
    pub fn new(a1: F, a2: F, a3: F, a4: F, a6: F) -> Self {
        Self { a1, a2, a3, a4, a6 }
    }

    /// Construct the **short Weierstrass** curve  `y² = x³ + ax + b`.
    ///
    /// Sets `a₁ = a₂ = a₃ = 0`,  `a₄ = a`,  `a₆ = b`.
    /// Only valid when `char(F) ≠ 2, 3`.
    pub fn new_short(a: F, b: F) -> Self {
        Self {
            a1: F::zero(),
            a2: F::zero(),
            a3: F::zero(),
            a4: a,
            a6: b,
        }
    }

    // -------------------------------------------------------------------
    // Invariants
    // -------------------------------------------------------------------

    /// Return the five a-invariants `[a₁, a₂, a₃, a₄, a₆]`.
    pub fn a_invariants(&self) -> [F; 5] {
        [
            self.a1.clone(),
            self.a2.clone(),
            self.a3.clone(),
            self.a4.clone(),
            self.a6.clone(),
        ]
    }

    /// Return the four b-invariants `[b₂, b₄, b₆, b₈]`.
    pub fn b_invariants(&self) -> [F; 4] {
        [self.b2(), self.b4(), self.b6(), self.b8()]
    }

    // -------------------------------------------------------------------
    // Curve predicates
    // -------------------------------------------------------------------

    /// Check whether the affine point `(x, y)` satisfies the curve equation.
    ///
    /// ```text
    /// y² + a₁xy + a₃y  ==  x³ + a₂x² + a₄x + a₆
    /// ```
    pub fn contains(&self, x: &F, y: &F) -> bool {
        let lhs = {
            let y2 = <F as FieldOps>::square(y);
            let a1xy = <F as FieldOps>::mul(&self.a1, &<F as FieldOps>::mul(x, y));
            let a3y = <F as FieldOps>::mul(&self.a3, y);
            <F as FieldOps>::add(&<F as FieldOps>::add(&y2, &a1xy), &a3y)
        };

        let rhs = {
            let x2 = <F as FieldOps>::square(x);
            let x3 = <F as FieldOps>::mul(&x2, x);
            let a2x2 = <F as FieldOps>::mul(&self.a2, &x2);
            let a4x = <F as FieldOps>::mul(&self.a4, x);
            <F as FieldOps>::add(
                &<F as FieldOps>::add(&<F as FieldOps>::add(&x3, &a2x2), &a4x),
                &self.a6,
            )
        };

        lhs == rhs
    }

    // -------------------------------------------------------------------
    // Discriminant helpers  (Silverman, §III.1)
    // -------------------------------------------------------------------

    /// `b₂ = a₁² + 4a₂`.
    pub fn b2(&self) -> F {
        let a1_sq = <F as FieldOps>::square(&self.a1);
        let four_a2 = <F as FieldOps>::double(&<F as FieldOps>::double(&self.a2));
        <F as FieldOps>::add(&a1_sq, &four_a2)
    }

    /// `b₄ = a₁a₃ + 2a₄`.
    pub fn b4(&self) -> F {
        let a1a3 = <F as FieldOps>::mul(&self.a1, &self.a3);
        let two_a4 = <F as FieldOps>::double(&self.a4);
        <F as FieldOps>::add(&a1a3, &two_a4)
    }

    /// `b₆ = a₃² + 4a₆`.
    pub fn b6(&self) -> F {
        let a3_sq = <F as FieldOps>::square(&self.a3);
        let four_a6 = <F as FieldOps>::double(&<F as FieldOps>::double(&self.a6));
        <F as FieldOps>::add(&a3_sq, &four_a6)
    }

    /// `b₈ = a₁²a₆ + 4a₂a₆ − a₁a₃a₄ + a₂a₃² − a₄²`.
    pub fn b8(&self) -> F {
        let a1sq_a6 = <F as FieldOps>::mul(&<F as FieldOps>::square(&self.a1), &self.a6);
        let four_a2_a6 = <F as FieldOps>::double(&<F as FieldOps>::double(&<F as FieldOps>::mul(
            &self.a2, &self.a6,
        )));
        let a1_a3_a4 = <F as FieldOps>::mul(&<F as FieldOps>::mul(&self.a1, &self.a3), &self.a4);
        let a2_a3sq = <F as FieldOps>::mul(&self.a2, &<F as FieldOps>::square(&self.a3));
        let a4sq = <F as FieldOps>::square(&self.a4);

        let sum = <F as FieldOps>::add(&a1sq_a6, &four_a2_a6);
        let sum = <F as FieldOps>::sub(&sum, &a1_a3_a4);
        let sum = <F as FieldOps>::add(&sum, &a2_a3sq);
        <F as FieldOps>::sub(&sum, &a4sq)
    }

    /// The discriminant  `Δ = −b₂²b₈ − 8b₄³ − 27b₆² + 9b₂b₄b₆`.
    ///
    /// The curve is non-singular if and only if `Δ ≠ 0`.
    /// (Only meaningful when `char(F) ≠ 2, 3`.)
    pub fn discriminant(&self) -> F {
        let b2 = self.b2();
        let b4 = self.b4();
        let b6 = self.b6();
        let b8 = self.b8();

        // Helper: build small integer constants from repeated doubling / adding
        let two = <F as FieldOps>::double(&F::one());
        let three = <F as FieldOps>::add(&two, &F::one());
        let eight = <F as FieldOps>::double(&<F as FieldOps>::double(&two));
        let nine = <F as FieldOps>::square(&three);
        let twentyseven = <F as FieldOps>::mul(&nine, &three);

        // −b₂²b₈
        let term1 =
            <F as FieldOps>::negate(&<F as FieldOps>::mul(&<F as FieldOps>::square(&b2), &b8));

        // −8b₄³
        let b4_cubed = <F as FieldOps>::mul(&<F as FieldOps>::square(&b4), &b4);
        let term2 = <F as FieldOps>::negate(&<F as FieldOps>::mul(&eight, &b4_cubed));

        // −27b₆²
        let term3 = <F as FieldOps>::negate(&<F as FieldOps>::mul(
            &twentyseven,
            &<F as FieldOps>::square(&b6),
        ));

        // 9b₂b₄b₆
        let term4 = <F as FieldOps>::mul(
            &nine,
            &<F as FieldOps>::mul(&b2, &<F as FieldOps>::mul(&b4, &b6)),
        );

        <F as FieldOps>::add(
            &<F as FieldOps>::add(&term1, &term2),
            &<F as FieldOps>::add(&term3, &term4),
        )
    }
}

impl<F: FieldOps> Curve for WeierstrassCurve<F> {
    type BaseField = F;
    type Point = AffinePoint<F>;

    fn is_on_curve(&self, point: &Self::Point) -> bool {
        if point.infinity {
            true
        } else {
            self.contains(&point.x, &point.y)
        }
    }

    fn random_point(&self) -> Self::Point {
        todo!()
    }

    fn j_invariant(&self) -> u64 {
        todo!()
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        WeierstrassCurve::a_invariants(self).to_vec()
    }
}
