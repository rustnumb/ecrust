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
        assert!(Self::is_smooth(&a1, &a2, &a3, &a4, &a6));

        Self {
            a1,
            a2,
            a3,
            a4,
            a6
        }
    }

    /// Construct the **short Weierstrass** curve  `y² = x³ + ax + b`.
    ///
    /// Sets `a₁ = a₂ = a₃ = 0`,  `a₄ = a`,  `a₆ = b`.
    /// Only valid when `char(F) ≠ 2, 3`.
    pub fn new_short(a: F, b: F) -> Self {
        assert!(F::characteristic()[0] > 3);
        assert!(Self::is_smooth(&F::zero(), &F::zero(), &F::zero(), &a, &b));

        Self {
            a1: F::zero(),
            a2: F::zero(),
            a3: F::zero(),
            a4: a,
            a6: b,
        }
    }

    pub fn is_smooth(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> bool {
        discriminant_from_coeffs::<F>(&a1, &a2, &a3, &a4, &a6) != F::zero()
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
            let a1xy = self.a1 * *x * *y;
            let a3y = self.a3 * *y;
            y2 + a1xy + a3y
        };

        let rhs = {
            let x2 = <F as FieldOps>::square(x);
            let x3 = *x * x2;
            let a2x2 = self.a2 * x2;
            let a4x = self.a4 * *x;
            x3 + a2x2 + a4x + self.a6
        };

        lhs == rhs
    }
}

// -------------------------------------------------------------------
// Discriminant private helpers  (Silverman, §III.1)
// -------------------------------------------------------------------

/// `b₂ = a₁² + 4a₂`.
fn b2_from_coeffs<F: FieldOps>(a1: &F, a2: &F) -> F {
    let a1_sq = <F as FieldOps>::square(&a1);
    let four_a2 = <F as FieldOps>::double(&<F as FieldOps>::double(&a2));
    a1_sq + four_a2
}

/// `b₄ = a₁a₃ + 2a₄`.
fn b4_from_coeffs<F: FieldOps>(a1: &F, a3: &F, a4: &F) -> F {
    let a1a3 = *a1 * *a3;
    let two_a4 = <F as FieldOps>::double(&a4);
    a1a3 + two_a4
}

/// `b₆ = a₃² + 4a₆`.
fn b6_from_coeffs<F: FieldOps>(a3: &F, a6: &F) -> F {
    let a3_sq = <F as FieldOps>::square(&a3);
    let four_a6 = <F as FieldOps>::double(&<F as FieldOps>::double(&a6));
    a3_sq + four_a6
}


/// `b₈ = a₁²a₆ + 4a₂a₆ − a₁a₃a₄ + a₂a₃² − a₄²`.
fn b8_from_coeffs<F: FieldOps>(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> F {
    let a1sq_a6 = <F as FieldOps>::square(&a1) * *a6;
    let four_a2_a6 = <F as FieldOps>::double(&<F as FieldOps>::double(&(*a2 * *a6) ));
    let a1_a3_a4 = *a1 * *a3 * *a4;
    let a2_a3sq = *a2 * <F as FieldOps>::square(&a3);
    let a4sq = <F as FieldOps>::square(&a4);

    a1sq_a6 + four_a2_a6 - a1_a3_a4 + a2_a3sq - a4sq
}

/// The discriminant  `Δ = −b₂²b₈ − 8b₄³ − 27b₆² + 9b₂b₄b₆`.
///
/// The curve is non-singular if and only if `Δ ≠ 0`.
/// (Only meaningful when `char(F) ≠ 2, 3`.)
fn discriminant_from_coeffs<F: FieldOps>(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> F {
    let b2 = b2_from_coeffs::<F>(&a1, &a2);
    let b4 = b4_from_coeffs::<F>(&a1, &a3, &a4);
    let b6 = b6_from_coeffs::<F>(&a3, &a6);
    let b8 = b8_from_coeffs::<F>(&a1, &a2, &a3, &a4, &a6);

    // Helper: build small integer constants from repeated doubling / adding
    let two = <F as FieldOps>::double(&F::one());
    let three = two + F::one();
    let eight = <F as FieldOps>::double(&<F as FieldOps>::double(&two));
    let nine = <F as FieldOps>::square(&three);
    let twentyseven = nine * three;

    // −b₂²b₈
    let term1 = - <F as FieldOps>::square(&b2) * b8;

    // −8b₄³
    let b4_cubed = <F as FieldOps>::square(&b4) * b4;
    let term2 = - eight * b4_cubed;

    // −27b₆²
    let term3 = - twentyseven * <F as FieldOps>::square(&b6);

    // 9b₂b₄b₆
    let term4 = nine * b2 * b4 * b6;

    term1 + term2 + term3 + term4
}

/// `c_4 = b_2² - 24b_4`.
fn c4_from_coeffs<F: FieldOps>(a1: &F, a2: &F, a3: &F, a4: &F) -> F {
    let b2 = b2_from_coeffs::<F>(&a1, &a2);
    let b4 = b4_from_coeffs::<F>(&a1, &a3, &a4);
    let b2_sq = <F as FieldOps>::square(&b2);

    // Helper: build small integer constants from repeated doubling / adding
    let two = <F as FieldOps>::double(&F::one());
    let three = two + F::one();
    let eight = <F as FieldOps>::double(&<F as FieldOps>::double(&two));
    let twentyfour = eight * three;

    b2_sq - twentyfour * b4
}

fn j_inv_from_coeffs<F: FieldOps>(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> F {
    let c4 = c4_from_coeffs::<F>(&a1, &a2, &a3, &a4);
    let c4_sq = <F as FieldOps>::square(&c4);
    let c4_cubed = c4 * c4_sq;

    let delta = discriminant_from_coeffs::<F>(&a1, &a2, &a3, &a4, &a6);
    let delta_inv = <F as FieldOps>::invert(&delta).unwrap();
    c4_cubed * delta_inv
}


// -------------------------------------------------------------------
// Invariants attached to the model
// -------------------------------------------------------------------

impl<F: FieldOps> WeierstrassCurve<F> {
    pub fn b2(&self) -> F {
        b2_from_coeffs::<F>(&self.a1, &self.a2)
    }
    pub fn b4(&self) -> F {
        b4_from_coeffs::<F>(&self.a1, &self.a3, &self.a4)
    }
    pub fn b6(&self) -> F {
        b6_from_coeffs::<F>(&self.a3, &self.a6)
    }
    pub fn b8(&self) -> F {
        b8_from_coeffs::<F>(&self.a1, &self.a2, &self.a3, &self.a4, &self.a6)    
    }
    pub fn discriminant(&self) -> F {
        discriminant_from_coeffs::<F>(&self.a1, &self.a2, &self.a3, &self.a4, &self.a6)
    }
}


// -------------------------------------------------------------------
// Curve predicates
// -------------------------------------------------------------------

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

    fn j_invariant(&self) -> F {
        j_inv_from_coeffs(&self.a1, &self.a2, &self.a3, &self.a4, &self.a6)
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        WeierstrassCurve::a_invariants(self).to_vec()
    }
}
