//! Elliptic curve definition in Jacobi intersection form.
//!
//! # Equation
//!
//! $$
//! s^2 + c^2 = 1,
//! \quad
//! a s^2 + d^2 = 1
//! $$
//!
//! over a field $F$ with $\mathrm{char}(F) \ne 2$.
//!
//! The EFD records that this model is birationally equivalent to the Weierstrass
//! curve
//!
//! ```text
//! y² = x³ + (2-a)x² + (1-a)x.
//! ```
//!

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::field_ops::{FieldOps, FieldRandom};
use fp::{ref_field_impl, ref_field_trait_impl};
use rand::{CryptoRng, Rng};

use crate::curve_ops::Curve;
use crate::curve_weierstrass::WeierstrassCurve;
use crate::point_jacobi_intersection::JacobiIntersectionPoint;

/// A Jacobi-intersection curve.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct JacobiIntersectionCurve<F: FieldOps> {
    /// The variable a in the definition of the curve
    pub a: F,
}

ref_field_impl! {
    impl<F> JacobiIntersectionCurve<F> {
        /// Construct a Jacobi intersection from its parameter `a`.
        pub fn new(a: F) -> Self {
            assert!(F::characteristic()[0] != 2, "Jacobi intersections require char(F) != 2");
            assert!(Self::is_smooth(&a), "singular Jacobi intersection");
            Self { a }
        }

        /// Smoothness requirement: `a != 0` and `a != 1`.
        pub fn is_smooth(a: &F) -> bool {
            let zero = F::zero();
            let one = F::one();
            a != &zero && a != &one
        }

        /// Check both defining quadrics.
        pub fn contains(&self, s: &F, c: &F, d: &F) -> bool {
            let s2 = <F as FieldOps>::square(s);
            let c2 = <F as FieldOps>::square(c);
            let d2 = <F as FieldOps>::square(d);

            let lhs1 = &s2 + &c2;
            let lhs2 = &(&self.a * &s2) + &d2;
            let one = F::one();

            lhs1 == one && lhs2 == one
        }

        pub fn a_invariants(&self) -> [F; 1] {
            [self.a.clone()]
        }

        /// Birationally equivalent Weierstrass model
        /// `y² = x³ + (2-a)x² + (1-a)x`.
        pub fn to_weierstrass_curve(&self) -> WeierstrassCurve<F> {
            let one = F::one();
            let two = <F as FieldOps>::double(&one);

            let a2 = &two - &self.a;
            let a4 = &one - &self.a;

            WeierstrassCurve::new(F::zero(), a2, F::zero(), a4, F::zero())
        }
    }
}

ref_field_impl!{
    impl<F: FieldOps + FieldRandom> JacobiIntersectionCurve<F> {
        /// Sample a random affine point on this Jacobi‑intersection curve using the
        /// provided RNG.
        ///
        /// The method repeatedly samples `s` and then solves the defining quadrics
        /// for `c` and `d` by square‑root extraction, returning a point
        /// `(s, c, d)` on the curve.
        pub fn random_point(&self, rng: &mut (impl CryptoRng + Rng)) -> JacobiIntersectionPoint<F> {
            loop {
                let s = F::random(rng);
                let s2 = <F as FieldOps>::square(&s);

                let c2 = &F::one() - &s2;
                let d2 = &F::one() - &(&self.a * &s2);

                if let (Some(c), Some(d)) = (c2.sqrt().into_option(), d2.sqrt().into_option()) {
                    let p = JacobiIntersectionPoint::new(s, c, d);
                    debug_assert!(self.is_on_curve(&p));
                    return p;
                }
            }
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + FieldRandom> Curve for JacobiIntersectionCurve<F> {
        type BaseField = F;
        type Point = JacobiIntersectionPoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            self.contains(&point.s, &point.c, &point.d)
        }

        fn random_point(&self, rng: &mut (impl CryptoRng + Rng)) -> Self::Point {
            self.random_point(rng)
        }

        fn j_invariant(&self) -> F {
            self.to_weierstrass_curve().j_invariant()
        }

        fn a_invariants(&self) -> Vec<Self::BaseField> {
            JacobiIntersectionCurve::a_invariants(self).to_vec()
        }
    }
}

// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for JacobiIntersectionCurve<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            a: F::conditional_select(&a.a, &b.a, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.a, &other.a, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.a, &mut b.a, choice);
    }
}

impl<F> ConstantTimeEq for JacobiIntersectionCurve<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.a.ct_eq(&other.a)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}