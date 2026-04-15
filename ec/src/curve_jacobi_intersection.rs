//! Elliptic curve definition in Jacobi intersection form.
//!
//! # Equation
//!
//! ```text
//! s² + c² = 1,
//! a s² + d² = 1
//! ```
//!
//! over a field `F` of characteristic different from `2`.
//!
//! The EFD records that this model is birationally equivalent to the Weierstrass
//! curve
//!
//! ```text
//! y² = x³ + (2-a)x² + (1-a)x.
//! ```
//!
use core::fmt;
use fp::field_ops::{FieldOps, FieldRandom};

use crate::curve_ops::Curve;
use crate::curve_weierstrass::WeierstrassCurve;
use crate::point_jacobi_intersection::JacobiIntersectionPoint;

/// A Jacobi-intersection curve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JacobiIntersectionCurve<F: FieldOps> {
    pub a: F,
}

impl<F> fmt::Display for JacobiIntersectionCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(
                f,
                "JacobiIntersectionCurve {{\n  s^2 + c^2 = 1\n  a s^2 + d^2 = 1\n  a = {}\n}}",
                self.a
            )
        } else {
            write!(f, "s^2 + c^2 = 1, {} s^2 + d^2 = 1", self.a)
        }
    }
}

impl<F: FieldOps + FieldRandom> JacobiIntersectionCurve<F> {
    /// Construct a Jacobi intersection from its parameter `a`.
    pub fn new(a: F) -> Self {
        assert!(F::characteristic()[0] != 2, "Jacobi intersections require char(F) != 2");
        assert!(Self::is_smooth(&a), "singular Jacobi intersection");
        Self { a }
    }

    /// Smoothness requirement: `a != 0` and `a != 1`.
    pub fn is_smooth(a: &F) -> bool {
        *a != F::zero() && *a != F::one()
    }

    /// Check both defining quadrics.
    pub fn contains(&self, s: &F, c: &F, d: &F) -> bool {
        let s2 = <F as FieldOps>::square(s);
        let c2 = <F as FieldOps>::square(c);
        let d2 = <F as FieldOps>::square(d);

        s2 + c2 == F::one() && self.a * s2 + d2 == F::one()
    }

    pub fn a_invariants(&self) -> [F; 1] {
        [self.a]
    }

    pub fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> JacobiIntersectionPoint<F> {
        loop {
            let s = F::random(rng);
            let s2 = <F as FieldOps>::square(&s);

            let c2 = F::one() - s2;
            let d2 = F::one() - self.a * s2;

            if let (Some(c), Some(d)) = (c2.sqrt().into_option(), d2.sqrt().into_option()) {
                let p = JacobiIntersectionPoint::new(s, c, d);
                debug_assert!(self.is_on_curve(&p));
                return p;
            }
        }
    }

    /// Birationally equivalent Weierstrass model
    /// $y^2 = x^3 + (2-a)x^2 + (1-a)x$.
    pub fn to_weierstrass_curve(&self) -> WeierstrassCurve<F> {
        let two = <F as FieldOps>::double(&F::one());
        WeierstrassCurve::new(F::zero(), two - self.a, F::zero(), F::one() - self.a, F::zero())
    }
}

impl<F: FieldOps + FieldRandom> Curve for JacobiIntersectionCurve<F> {
    type BaseField = F;
    type Point = JacobiIntersectionPoint<F>;

    fn is_on_curve(&self, point: &Self::Point) -> bool {
        self.contains(&point.s, &point.c, &point.d)
    }

    fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
        JacobiIntersectionCurve::random_point(self, rng)
    }

    fn j_invariant(&self) -> F {
        self.to_weierstrass_curve().j_invariant()
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        JacobiIntersectionCurve::a_invariants(self).to_vec()
    }
}
