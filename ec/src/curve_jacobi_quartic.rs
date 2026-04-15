//! Elliptic curve definition in (extended) Jacobi quartic form.
//!
//! # Equation
//!
//! ```text
//! y² = d x⁴ + 2 a x² + 1
//! ```
//!
//! over a field `F` of characteristic different from `2`.
//!
//! This follows the model studied in Hisil–Wong–Carter–Dawson,
//! *Jacobi Quartic Curves Revisited* (2009), which treats the more general
//! “extended Jacobi quartic” family with arbitrary `a` and `d` satisfying
//! `d(a²-d) ≠ 0`.

use fp::field_ops::FieldOps;

use crate::curve_ops::Curve;
use crate::point_jacobi_quartic::JacobiQuarticPoint;

/// A Jacobi quartic curve
///
/// ```text
/// y² = d x⁴ + 2 a x² + 1
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JacobiQuarticCurve<F: FieldOps> {
    pub a: F,
    pub d: F,
}

impl<F: FieldOps> JacobiQuarticCurve<F> {
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
        a2 != *d
    }

    /// Check whether `(x, y)` lies on `y² = d x⁴ + 2 a x² + 1`.
    pub fn contains(&self, x: &F, y: &F) -> bool {
        let x2 = <F as FieldOps>::square(x);
        let x4 = <F as FieldOps>::square(&x2);
        let y2 = <F as FieldOps>::square(y);
        let two = <F as FieldOps>::double(&F::one());

        y2 == self.d * x4 + two * self.a * x2 + F::one()
    }

    /// Return `[a, d]`.
    pub fn a_invariants(&self) -> [F; 2] {
        [self.a, self.d]
    }

    /// Return the affine identity `(0, 1)`.
    pub fn neutral_point(&self) -> JacobiQuarticPoint<F> {
        JacobiQuarticPoint::identity()
    }
}

impl<F: FieldOps> Curve for JacobiQuarticCurve<F> {
    type BaseField = F;
    type Point = JacobiQuarticPoint<F>;

    fn is_on_curve(&self, point: &Self::Point) -> bool {
        self.contains(&point.x, &point.y)
    }

    fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
        todo!()
    }

    /// The paper gives
    ///
    /// ```text
    /// j = 64 d^{-1} (a²-d)^{-2} (a²+3d)³.
    /// ```
    ///
    fn j_invariant(&self) -> F {
        let a2 = <F as FieldOps>::square(&self.a);
        let three = <F as FieldOps>::double(&F::one()) + F::one();

        let eight = <F as FieldOps>::double(
            &<F as FieldOps>::double(
                &<F as FieldOps>::double(&F::one())
            )
        );
        let sixty_four = eight * eight; // 8 * 8 = 64

        let num_base = a2 + three * self.d;
        let num = sixty_four * num_base * <F as FieldOps>::square(&num_base);

        let diff = a2 - self.d;
        let denom = self.d * <F as FieldOps>::square(&diff);
        let denom_inv = denom
            .invert()
            .into_option()
            .expect("Jacobi quartic j-invariant denominator must be invertible");

        num * denom_inv
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        JacobiQuarticCurve::a_invariants(self).to_vec()
    }
}
