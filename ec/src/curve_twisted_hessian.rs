//! Elliptic curve definition in twisted Hessian form.
//!
//! # Equation
//!
//! We use the twisted Hessian model
//!
//! $$
//! aX^3 + Y^3 + Z^3 = 3dXYZ,
//! $$
//!
//! together with the affine chart
//!
//! $$
//! ax^3 + y^3 + 1 = 3dxy.
//! $$
//!
//! The neutral element is
//!
//! $$O = (0 : -1 : 1).$$
//!
//! # Smoothness
//!
//! The twisted Hessian curve is nonsingular when
//!
//! $$a \neq 0 \quad\text{and}\quad d^3 \neq a.$$
//!
//! # References
//!
//! - Thomas Decru and Sabrina Kunzweiler,
//!   *Tripling on Hessian curves via isogeny decomposition*, 2026.
//! - Daniel J. Bernstein, Chitchanok Chuengsatiansup,
//!   David Kohel, and Tanja Lange,
//!   *Twisted Hessian curves*, LATINCRYPT 2015.

use core::fmt;

use fp::field_ops::{FieldOps, FieldRandom};
use fp::{ref_field_impl, ref_field_trait_impl};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_ops::Curve;
use crate::point_twisted_hessian::TwistedHessianPoint;

/// A twisted Hessian curve
///
/// $$aX^3 + Y^3 + Z^3 = 3dXYZ.$$
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct TwistedHessianCurve<F: FieldOps> {
    /// The twisting parameter `a`.
    pub a: F,
    /// The Hessian parameter `d`.
    pub d: F,
}

impl<F> fmt::Display for TwistedHessianCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(
                f,
                "TwistedHessianCurve {{\n  aX^3 + Y^3 + Z^3 = 3dXYZ\n  a = {}\n  d = {}\n}}",
                self.a, self.d
            )
        } else {
            write!(f, "({})X^3 + Y^3 + Z^3 = 3({})XYZ", self.a, self.d)
        }
    }
}

ref_field_impl! {
    impl<F: FieldOps + FieldRandom> TwistedHessianCurve<F> {
        /// Construct a twisted Hessian curve.
        pub fn new(a: F, d: F) -> Self {
            assert!(Self::is_smooth(&a, &d), "singular twisted Hessian curve");
            Self { a, d }
        }

        /// Construct a twisted Hessian curve in twisted normal form (`d = 1`).
        pub fn new_normal_form(a: F) -> Self {
            Self::new(a, F::one())
        }

        /// Return `true` if the twisted Hessian discriminant is nonzero.
        pub fn is_smooth(a: &F, d: &F) -> bool {
            if bool::from(a.is_zero()) {
                return false;
            }

            let d2 = <F as FieldOps>::square(d);
            let d3 = d * &d2;
            d3 != a.clone()
        }

        /// Check whether the affine point `(x, y)` satisfies
        ///
        /// $$ax^3 + y^3 + 1 = 3dxy.$$
        pub fn contains_affine(&self, x: &F, y: &F) -> bool {
            let x2 = <F as FieldOps>::square(x);
            let y2 = <F as FieldOps>::square(y);
            let x3 = x * &x2;
            let y3 = y * &y2;

            let lhs = &(&self.a * &x3) + &(&y3 + &F::one());
            let rhs = &(&F::from_u64(3) * &self.d) * &(x * y);
            lhs == rhs
        }

        /// Check whether the projective point `(X:Y:Z)` satisfies
        ///
        /// $$aX^3 + Y^3 + Z^3 = 3dXYZ.$$
        pub fn contains_projective(&self, x: &F, y: &F, z: &F) -> bool {
            let x2 = <F as FieldOps>::square(x);
            let y2 = <F as FieldOps>::square(y);
            let z2 = <F as FieldOps>::square(z);
            let x3 = x * &x2;
            let y3 = y * &y2;
            let z3 = z * &z2;

            let lhs = &(&self.a * &x3) + &(&y3 + &z3);
            let rhs = &(&F::from_u64(3) * &self.d) * &(x * &(y * z));
            lhs == rhs
        }

        /// Return `[a, d]`.
        pub fn a_invariants(&self) -> [F; 2] {
            [self.a.clone(), self.d.clone()]
        }

        /// Return the neutral element `(0:-1:1)`.
        pub fn neutral_point(&self) -> TwistedHessianPoint<F> {
            TwistedHessianPoint::identity()
        }

        /// Best-effort random point sampling in the affine chart `Z = 1`.
        pub fn random_point(
            &self,
            rng: &mut (impl rand::CryptoRng + rand::Rng),
        ) -> TwistedHessianPoint<F> {
            loop {
                let x = F::random(rng);
                let y = F::random(rng);
                if self.contains_affine(&x, &y) {
                    let p = TwistedHessianPoint::from_affine(x, y);
                    debug_assert!(self.is_on_curve(&p));
                    return p;
                }
            }
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + FieldRandom> Curve for TwistedHessianCurve<F> {
        type BaseField = F;
        type Point = TwistedHessianPoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            self.contains_projective(&point.x, &point.y, &point.z)
        }

        fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
            TwistedHessianCurve::random_point(self, rng)
        }

        fn j_invariant(&self) -> F {
            // Corollary 2.10 in Decru--Kunzweiler (2026):
            // j(H_{a,d}) = a^{-1} * [ 3 d (8a + d^3) / (d^3 - a) ]^3.
            let d2 = <F as FieldOps>::square(&self.d);
            let d3 = &self.d * &d2;

            let eight_a = &F::from_u64(8) * &self.a;
            let inner_num = &(&F::from_u64(3) * &self.d) * &(&eight_a + &d3);
            let inner_den = &d3 - &self.a;
            let inner = &inner_num
                * &inner_den
                    .invert()
                    .into_option()
                    .expect("twisted Hessian j-invariant denominator must be invertible");
            let inner_cubed = &inner * &<F as FieldOps>::square(&inner);

            &self.a
                .invert()
                .into_option()
                .expect("a must be invertible on a nonsingular twisted Hessian curve")
                * &inner_cubed
        }

        fn a_invariants(&self) -> Vec<Self::BaseField> {
            TwistedHessianCurve::a_invariants(self).to_vec()
        }
    }
}

impl<F> ConditionallySelectable for TwistedHessianCurve<F>
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
        self.a.conditional_assign(&other.a, choice);
        self.d.conditional_assign(&other.d, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.a, &mut b.a, choice);
        F::conditional_swap(&mut a.d, &mut b.d, choice);
    }
}

impl<F> ConstantTimeEq for TwistedHessianCurve<F>
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
