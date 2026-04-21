//! Elliptic curve definition (general Weierstrass form).
//!
//! # Equation
//!
//! The **general Weierstrass** equation over a field $F$ is
//!
//! $$
//! y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6
//! $$
//!
//! This form is valid in *any* characteristic, including characteristic $2$.
//!
//! # Short Weierstrass specialisation
//!
//! When $\mathrm{char}(F) \ne 2, 3$ the curve can be brought to the simpler
//!
//! $$
//! y^2 = x^3 + ax + b \quad (a_1 = a_2 = a_3 = 0,\; a_4 = a,\; a_6 = b)
//! $$
//!
//! via the convenience constructor [`WeierstrassCurve::new_short`].

use core::fmt;
use fp::field_ops::{FieldOps, FieldRandom};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::{ref_field_fns, ref_field_impl, ref_field_trait_impl};

use crate::curve_ops::Curve;
use crate::point_weierstrass::AffinePoint;

/// An elliptic curve in general Weierstrass form over a field $F$.
///
/// $$
/// y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6
/// $$
///
/// All five coefficients are stored explicitly; the short Weierstrass
/// case simply has $a_1 = a_2 = a_3 = 0$.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct WeierstrassCurve<F: FieldOps> {
    /// a-invariant
    pub a1: F,
    /// a-invariant
    pub a2: F,
    /// a-invariant
    pub a3: F,
    /// a-invariant
    pub a4: F,
    /// a-invariant
    pub a6: F,
}

impl<F> fmt::Display for WeierstrassCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(
                f,
                "WeierstrassCurve {{\n  y^2 + ({})xy + ({})y = x^3 + ({})x^2 + ({})x + ({})\n}}",
                self.a1, self.a3, self.a2, self.a4, self.a6
            )
        } else {
            write!(
                f,
                "y^2 + ({})xy + ({})y = x^3 + ({})x^2 + ({})x + ({})",
                self.a1, self.a3, self.a2, self.a4, self.a6
            )
        }
    }
}

ref_field_impl! {
    impl<F> WeierstrassCurve<F> {
        // -------------------------------------------------------------------
        // Constructors
        // -------------------------------------------------------------------

        /// Constructs a curve from the five general Weierstrass coefficients
        /// $(a_1, a_2, a_3, a_4, a_6)$.
        ///
        /// The curve is required to be smooth, i.e. its discriminant $\Delta \ne 0$.
        ///
        /// # Panics
        ///
        /// Panics if the curve is singular.
        pub fn new(a1: F, a2: F, a3: F, a4: F, a6: F) -> Self {
            assert!(Self::is_smooth(&a1, &a2, &a3, &a4, &a6));
            Self { a1, a2, a3, a4, a6 }
        }

        /// Constructs the **short Weierstrass** curve
        ///
        /// $$
        /// y^2 = x^3 + ax + b.
        /// $$
        ///
        /// This sets $a_1 = a_2 = a_3 = 0$, $a_4 = a$, and $a_6 = b$.
        ///
        /// This form is only valid when $\mathrm{char}(F) \ne 2, 3$.
        ///
        /// # Panics
        ///
        /// Panics if the characteristic is $\le 3$ or if the curve is singular.
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

        /// Returns `true` if and only if the curve defined by the given
        /// $a$-invariants is smooth.
        ///
        /// This is equivalent to checking that the discriminant satisfies
        /// $\Delta \ne 0$.
        pub fn is_smooth(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> bool {
            discriminant_from_coeffs::<F>(a1, a2, a3, a4, a6) != F::zero()
        }

        // -------------------------------------------------------------------
        // Invariants
        // -------------------------------------------------------------------

        /// Returns the five $a$-invariants $[a_1, a_2, a_3, a_4, a_6]$.
        pub fn a_invariants(&self) -> [F; 5] {
            [
                self.a1.clone(),
                self.a2.clone(),
                self.a3.clone(),
                self.a4.clone(),
                self.a6.clone(),
            ]
        }

        /// Returns the four $b$-invariants $[b_2, b_4, b_6, b_8]$.
        pub fn b_invariants(&self) -> [F; 4] {
            [self.b2(), self.b4(), self.b6(), self.b8()]
        }

        // -------------------------------------------------------------------
        // Curve predicates
        // -------------------------------------------------------------------

        /// Checks whether the affine point $(x, y)$ satisfies the curve equation.
        ///
        /// $$
        /// y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6.
        /// $$
        pub fn contains(&self, x: &F, y: &F) -> bool {
            let lhs = {
                let y2 = <F as FieldOps>::square(y);
                let xy = x * y;
                let a1xy = &self.a1 * &xy;
                let a3y = &self.a3 * y;
                let tmp = &y2 + &a1xy;
                &tmp + &a3y
            };

            let rhs = {
                let x2 = <F as FieldOps>::square(x);
                let x3 = x * &x2;
                let a2x2 = &self.a2 * &x2;
                let a4x = &self.a4 * x;
                let tmp1 = &x3 + &a2x2;
                let tmp2 = &tmp1 + &a4x;
                &tmp2 + &self.a6
            };

            lhs == rhs
        }
    }
}

ref_field_impl!{
    impl<F: FieldOps + FieldRandom> WeierstrassCurve<F> {
        /// Sample a random affine point on this curve using the provided RNG.
        ///
        /// This currently uses a square-root-based construction and is implemented
        /// only for odd characteristic. It returns a finite affine point `P` such
        /// that `self.is_on_curve(&P)` holds.
        pub fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> AffinePoint<F> {
            assert!(
                F::characteristic()[0] != 2,
                "random_point currently implemented for odd characteristic"
            );

            let one = F::one();
            let two = <F as FieldOps>::double(&one);
            let two_inv = <F as FieldOps>::invert(&two)
                .into_option()
                .expect("2 must be invertible");

            loop {
                let x = <F as FieldRandom>::random(rng);

                let a1x = &self.a1 * &x;
                let b = &a1x + &self.a3;

                let x2 = <F as FieldOps>::square(&x);
                let x3 = &x * &x2;
                let a2x2 = &self.a2 * &x2;
                let a4x = &self.a4 * &x;

                let rhs_tmp1 = &x3 + &a2x2;
                let rhs_tmp2 = &rhs_tmp1 + &a4x;
                let rhs = &rhs_tmp2 + &self.a6;

                let b_sq = <F as FieldOps>::square(&b);
                let two_rhs = <F as FieldOps>::double(&rhs);
                let four_rhs = <F as FieldOps>::double(&two_rhs);
                let disc = &b_sq + &four_rhs;

                if let Some(sqrt_disc) = disc.sqrt().into_option() {
                    let neg_b = -&b;
                    let sum = &neg_b + &sqrt_disc;
                    let y = &sum * &two_inv;
                    let p = AffinePoint::new(x, y);
                    debug_assert!(self.is_on_curve(&p));
                    return p;
                }
            }
        }
    }
}


// -------------------------------------------------------------------
// Discriminant private helpers  (Silverman, §III.1)
// -------------------------------------------------------------------

ref_field_fns! {
    /// `b₂ = a₁² + 4a₂`.
    fn b2_from_coeffs<F>(a1: &F, a2: &F) -> F {
        let a1_sq = <F as FieldOps>::square(a1);
        let two_a2 = <F as FieldOps>::double(a2);
        let four_a2 = <F as FieldOps>::double(&two_a2);
        &a1_sq + &four_a2
    }

    /// `b₄ = a₁a₃ + 2a₄`.
    fn b4_from_coeffs<F>(a1: &F, a3: &F, a4: &F) -> F {
        let a1a3 = a1 * a3;
        let two_a4 = <F as FieldOps>::double(a4);
        &a1a3 + &two_a4
    }

    /// `b₆ = a₃² + 4a₆`.
    fn b6_from_coeffs<F>(a3: &F, a6: &F) -> F {
        let a3_sq = <F as FieldOps>::square(a3);
        let two_a6 = <F as FieldOps>::double(a6);
        let four_a6 = <F as FieldOps>::double(&two_a6);
        &a3_sq + &four_a6
    }

    /// `b₈ = a₁²a₆ + 4a₂a₆ − a₁a₃a₄ + a₂a₃² − a₄²`.
    fn b8_from_coeffs<F>(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> F {
        let a1_sq = <F as FieldOps>::square(a1);
        let a1sq_a6 = &a1_sq * a6;

        let a2a6 = a2 * a6;
        let two_a2_a6 = <F as FieldOps>::double(&a2a6);
        let four_a2_a6 = <F as FieldOps>::double(&two_a2_a6);

        let a1a3 = a1 * a3;
        let a1_a3_a4 = &a1a3 * a4;

        let a3_sq = <F as FieldOps>::square(a3);
        let a2_a3sq = a2 * &a3_sq;

        let a4_sq = <F as FieldOps>::square(a4);

        let tmp1 = &a1sq_a6 + &four_a2_a6;
        let tmp2 = &tmp1 - &a1_a3_a4;
        let tmp3 = &tmp2 + &a2_a3sq;
        &tmp3 - &a4_sq
    }

    /// The discriminant  `Δ = −b₂²b₈ − 8b₄³ − 27b₆² + 9b₂b₄b₆`.
    ///
    /// The curve is non-singular if and only if `Δ ≠ 0`.
    /// (Only meaningful when `char(F) ≠ 2, 3`.)
    fn discriminant_from_coeffs<F>(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> F {
        let b2 = b2_from_coeffs::<F>(a1, a2);
        let b4 = b4_from_coeffs::<F>(a1, a3, a4);
        let b6 = b6_from_coeffs::<F>(a3, a6);
        let b8 = b8_from_coeffs::<F>(a1, a2, a3, a4, a6);

        let one = F::one();
        let two = <F as FieldOps>::double(&one);
        let three = &two + &one;
        let four = <F as FieldOps>::double(&two);
        let eight = <F as FieldOps>::double(&four);
        let nine = <F as FieldOps>::square(&three);
        let twentyseven = &nine * &three;

        let b2_sq = <F as FieldOps>::square(&b2);
        let neg_b2_sq = -&b2_sq;
        let term1 = &neg_b2_sq * &b8;

        let b4_sq = <F as FieldOps>::square(&b4);
        let b4_cubed = &b4_sq * &b4;
        let neg_eight = -&eight;
        let term2 = &neg_eight * &b4_cubed;

        let b6_sq = <F as FieldOps>::square(&b6);
        let neg_twentyseven = -&twentyseven;
        let term3 = &neg_twentyseven * &b6_sq;

        let b4b6 = &b4 * &b6;
        let b2b4b6 = &b2 * &b4b6;
        let term4 = &nine * &b2b4b6;

        let tmp1 = &term1 + &term2;
        let tmp2 = &tmp1 + &term3;
        &tmp2 + &term4
    }

    /// `c₄ = b₂² - 24b₄`.
    fn c4_from_coeffs<F>(a1: &F, a2: &F, a3: &F, a4: &F) -> F {
        let b2 = b2_from_coeffs::<F>(a1, a2);
        let b4 = b4_from_coeffs::<F>(a1, a3, a4);
        let b2_sq = <F as FieldOps>::square(&b2);

        let one = F::one();
        let two = <F as FieldOps>::double(&one);
        let three = &two + &one;
        let four = <F as FieldOps>::double(&two);
        let eight = <F as FieldOps>::double(&four);
        let twentyfour = &eight * &three;

        let twentyfour_b4 = &twentyfour * &b4;
        &b2_sq - &twentyfour_b4
    }

    fn j_inv_from_coeffs<F>(a1: &F, a2: &F, a3: &F, a4: &F, a6: &F) -> F {
        let c4 = c4_from_coeffs::<F>(a1, a2, a3, a4);
        let c4_sq = <F as FieldOps>::square(&c4);
        let c4_cubed = &c4 * &c4_sq;

        let delta = discriminant_from_coeffs::<F>(a1, a2, a3, a4, a6);
        let delta_inv = <F as FieldOps>::invert(&delta).unwrap();
        &c4_cubed * &delta_inv
    }
}

// -------------------------------------------------------------------
// Invariants attached to the model
// -------------------------------------------------------------------

ref_field_impl!{
    impl<F: FieldOps> WeierstrassCurve<F> {
        /// Returns the invariant $b_2 = a_1^2 + 4a_2$.
        pub fn b2(&self) -> F {
            b2_from_coeffs::<F>(&self.a1, &self.a2)
        }

        /// Returns the invariant $b_4 = a_1 a_3 + 2a_4$.
        pub fn b4(&self) -> F {
            b4_from_coeffs::<F>(&self.a1, &self.a3, &self.a4)
        }

        /// Returns the invariant $b_6 = a_3^2 + 4a_6$.
        pub fn b6(&self) -> F {
            b6_from_coeffs::<F>(&self.a3, &self.a6)
        }

        /// Returns the invariant
        ///
        /// $$
        /// b_8 = a_1^2 a_6 + 4 a_2 a_6 - a_1 a_3 a_4 + a_2 a_3^2 - a_4^2.
        /// $$
        pub fn b8(&self) -> F {
            b8_from_coeffs::<F>(&self.a1, &self.a2, &self.a3, &self.a4, &self.a6)
        }

        /// Returns the discriminant $\Delta$ of the curve.
        ///
        /// The curve is non-singular if and only if $\Delta \ne 0$.
        pub fn discriminant(&self) -> F {
            discriminant_from_coeffs::<F>(&self.a1, &self.a2, &self.a3, &self.a4, &self.a6)
        }
    }
}


// -------------------------------------------------------------------
// Curve predicates
// -------------------------------------------------------------------

ref_field_trait_impl!{
    impl<F: FieldOps + FieldRandom> Curve for WeierstrassCurve<F> {
    type BaseField = F;
    type Point = AffinePoint<F>;

    /// Returns `true` if the given point lies on the curve.
    ///
    /// The point at infinity is always considered to be on the curve.
    fn is_on_curve(&self, point: &Self::Point) -> bool {
        if point.infinity {
            true
        } else {
            self.contains(&point.x, &point.y)
        }
    }

    fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
        WeierstrassCurve::random_point(self, rng)
    }

    /// Returns the $j$-invariant of the curve.
    ///
    /// This is a complete invariant of elliptic curves over algebraically
    /// closed fields up to isomorphism.
    fn j_invariant(&self) -> F {
        j_inv_from_coeffs(&self.a1, &self.a2, &self.a3, &self.a4, &self.a6)
    }

    /// Returns the $a$-invariants as a vector.
    fn a_invariants(&self) -> Vec<Self::BaseField> {
        WeierstrassCurve::a_invariants(self).to_vec()
    }
}
}


// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for WeierstrassCurve<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            a1: F::conditional_select(&a.a1, &b.a1, choice),
            a2: F::conditional_select(&a.a2, &b.a2, choice),
            a3: F::conditional_select(&a.a3, &b.a3, choice),
            a4: F::conditional_select(&a.a4, &b.a4, choice),
            a6: F::conditional_select(&a.a6, &b.a6, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.a1, &other.a1, choice);
        F::conditional_assign(&mut self.a2, &other.a2, choice);
        F::conditional_assign(&mut self.a3, &other.a3, choice);
        F::conditional_assign(&mut self.a4, &other.a4, choice);
        F::conditional_assign(&mut self.a6, &other.a6, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.a1, &mut b.a1, choice);
        F::conditional_swap(&mut a.a2, &mut b.a2, choice);
        F::conditional_swap(&mut a.a3, &mut b.a3, choice);
        F::conditional_swap(&mut a.a4, &mut b.a4, choice);
        F::conditional_swap(&mut a.a6, &mut b.a6, choice);
    }
}

impl<F> ConstantTimeEq for WeierstrassCurve<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.a1.ct_eq(&other.a1)
            & self.a2.ct_eq(&other.a2)
            & self.a3.ct_eq(&other.a3)
            & self.a4.ct_eq(&other.a4)
            & self.a6.ct_eq(&other.a6)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}