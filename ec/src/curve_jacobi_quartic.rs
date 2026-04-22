//! Elliptic curve definition in (extended) Jacobi quartic form.
//!
//! # Equation
//!
//! The curve is given by
//!
//! $$
//! y^2 = d x^4 + 2 a x^2 + 1
//! $$
//!
//! over a field $F$ with $\mathrm{char}(F) \ne 2$.
//!
//! This follows the model studied in Hisil–Wong–Carter–Dawson,
//! *Jacobi Quartic Curves Revisited* (2009), which treats the more general
//! “extended Jacobi quartic” family with arbitrary `a` and `d` satisfying
//! `d(a²-d) ≠ 0`.

use core::fmt;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::field_ops::{FieldOps, FieldRandom};
use fp::{ref_field_impl, ref_field_trait_impl};

use crate::curve_ops::Curve;
use crate::point_jacobi_quartic::JacobiQuarticPoint;

/// A Jacobi quartic curve
///
/// ```text
/// y² = d x⁴ + 2 a x² + 1
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct JacobiQuarticCurve<F: FieldOps> {
    /// Invariant a in the definition
    pub a: F,
    /// Invariant d in the definition
    pub d: F,
}

impl<F> fmt::Display for JacobiQuarticCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(
                f,
                "JacobiQuarticCurve {{\n  y^2 = d x^4 + 2 a x^2 + 1\n  a = {}\n  d = {}\n}}",
                self.a, self.d
            )
        } else {
            write!(
                f,
                "y^2 = ({})x^4 + 2({})x^2 + 1",
                self.d, self.a
            )
        }
    }
}

ref_field_impl! {
    impl<F> JacobiQuarticCurve<F> {
        /// Construct a Jacobi quartic curve from `(a, d)`.
        pub fn new(a: F, d: F) -> Self {
            assert!(F::characteristic()[0] != 2, "Jacobi quartics require char(F) != 2");
            assert!(Self::is_smooth(&a, &d), "singular Jacobi quartic");
            Self { a, d }
        }

        /// Smoothness criterion from the discriminant
        /// `Δ = 256 d (a² - d)² ≠ 0`.
        pub fn is_smooth(a: &F, d: &F) -> bool {
            if bool::from(d.is_zero()) {
                return false;
            }

            let a2 = <F as FieldOps>::square(a);
            a2 != d.clone()
        }

        /// Check whether `(x, y)` lies on `y² = d x⁴ + 2 a x² + 1`.
        pub fn contains(&self, x: &F, y: &F) -> bool {
            let x2 = <F as FieldOps>::square(x);
            let x4 = <F as FieldOps>::square(&x2);
            let y2 = <F as FieldOps>::square(y);

            let one = F::one();
            let two = <F as FieldOps>::double(&one);

            let dx4 = &self.d * &x4;
            let ax2 = &self.a * &x2;
            let two_ax2 = &two * &ax2;

            let rhs_tmp = &dx4 + &two_ax2;
            let rhs = &rhs_tmp + &one;

            y2 == rhs
        }

        /// Return `[a, d]`.
        pub fn a_invariants(&self) -> [F; 2] {
            [self.a.clone(), self.d.clone()]
        }

        /// Return the affine identity `(0, 1)`.
        pub fn neutral_point(&self) -> JacobiQuarticPoint<F> {
            JacobiQuarticPoint::identity()
        }
    }
}

ref_field_impl!{
    impl<F: FieldOps + FieldRandom> JacobiQuarticCurve<F> {
        /// Sample a random affine point on this Jacobi quartic using the provided RNG.
        ///
        /// The method repeatedly samples `x`, evaluates the right-hand side of the
        /// quartic equation, and returns `(x, y)` when that value is a square in the
        /// base field.
        pub fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> JacobiQuarticPoint<F> {
            loop {
                let x = F::random(rng);
                let x2 = <F as FieldOps>::square(&x);
                let x4 = <F as FieldOps>::square(&x2);
                let rhs = &(&(&self.d * &x4)
                    + &(&<F as FieldOps>::double(&F::one()) * &(&self.a * &x2)))
                    + &F::one();

                if let Some(y) = rhs.sqrt().into_option() {
                    let p = JacobiQuarticPoint::new(x, y);
                    debug_assert!(self.is_on_curve(&p));
                    return p;
                }
            }
        }
    }
}


ref_field_trait_impl! {
    impl<F: FieldOps + FieldRandom> Curve for JacobiQuarticCurve<F> {
        type BaseField = F;
        type Point = JacobiQuarticPoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            self.contains(&point.x, &point.y)
        }

        fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
            self.random_point(rng)
        }

        /// The paper gives
        ///
        /// ```text
        /// j = 64 d^{-1} (a²-d)^{-2} (a²+3d)³.
        /// ```
        fn j_invariant(&self) -> F {
            let a2 = <F as FieldOps>::square(&self.a);

            let one = F::one();
            let two = <F as FieldOps>::double(&one);
            let three = &two + &one;

            let four = <F as FieldOps>::double(&two);
            let eight = <F as FieldOps>::double(&four);
            let sixty_four = &eight * &eight;

            let three_d = &three * &self.d;
            let num_base = &a2 + &three_d;
            let num_base_sq = <F as FieldOps>::square(&num_base);
            let num_base_cubed = &num_base * &num_base_sq;
            let num = &sixty_four * &num_base_cubed;

            let diff = &a2 - &self.d;
            let diff_sq = <F as FieldOps>::square(&diff);
            let denom = &self.d * &diff_sq;
            let denom_inv = <F as FieldOps>::invert(&denom)
                .into_option()
                .expect("Jacobi quartic j-invariant denominator must be invertible");

            &num * &denom_inv
        }

        fn a_invariants(&self) -> Vec<Self::BaseField> {
            JacobiQuarticCurve::a_invariants(self).to_vec()
        }
    }
}

// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for JacobiQuarticCurve<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            a: F::conditional_select(&a.a, &b.a, choice),
            d: F::conditional_select(&a.d, &b.d, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.a, &other.a, choice);
        F::conditional_assign(&mut self.d, &other.d, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.a, &mut b.a, choice);
        F::conditional_swap(&mut a.d, &mut b.d, choice);
    }
}

impl<F> ConstantTimeEq for JacobiQuarticCurve<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.a.ct_eq(&other.a) & self.d.ct_eq(&other.d)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}