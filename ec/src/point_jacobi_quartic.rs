//! Affine points on a Jacobi quartic curve.
//!
//! We use the affine model
//!
//! ```text
//! y² = d x⁴ + 2 a x² + 1
//! ```
//!
//! with neutral element `(0, 1)` and negation `-(x, y) = (-x, y)`. The
//! doubling formulas are taken from equations (9) and (10) in
//! *Jacobi Quartic Curves Revisited*; the general addition uses the affine
//! formulas (1) and (2).
//!
//! Important: these are affine formulas. Like the existing Edwards code in this
//! crate, they assume the denominators are invertible. For exceptional inputs
//! where the result leaves this affine chart, the code panics instead of trying
//! to model the desingularized points at infinity.

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_jacobi_quartic::JacobiQuarticCurve;
use crate::point_ops::{PointAdd, PointOps};
use fp::field_ops::FieldOps;

/// An affine point `(x, y)` on a Jacobi quartic curve.
#[derive(Debug, Clone, Copy)]
pub struct JacobiQuarticPoint<F: FieldOps> {
    pub x: F,
    pub y: F,
}

impl<F: FieldOps> PartialEq for JacobiQuarticPoint<F> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<F: FieldOps> Eq for JacobiQuarticPoint<F> {}

impl<F: FieldOps> JacobiQuarticPoint<F> {
    pub fn new(x: F, y: F) -> Self {
        Self { x, y }
    }

    /// The neutral element `(0, 1)`.
    pub fn identity() -> Self {
        Self { x: F::zero(), y: F::one() }
    }

    pub fn is_identity(&self) -> bool {
        self.x == F::zero() && self.y == F::one()
    }

    /// The affine order-2 point `(0, -1)`.
    pub fn order_two_point() -> Self {
        Self { x: F::zero(), y: -F::one() }
    }
}

impl<F> ConditionallySelectable for JacobiQuarticPoint<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            x: F::conditional_select(&a.x, &b.x, choice),
            y: F::conditional_select(&a.y, &b.y, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.x.conditional_assign(&other.x, choice);
        self.y.conditional_assign(&other.y, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.x, &mut b.x, choice);
        F::conditional_swap(&mut a.y, &mut b.y, choice);
    }
}

impl<F> ConstantTimeEq for JacobiQuarticPoint<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.x.ct_eq(&other.x) & self.y.ct_eq(&other.y)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}

impl<F: FieldOps> JacobiQuarticPoint<F> {
    /// Negation on a Jacobi quartic: `-(x, y) = (-x, y)`.
    pub fn negate(&self, _curve: &JacobiQuarticCurve<F>) -> Self {
        Self::new(-self.x, self.y)
    }

    /// Dedicated affine doubling from equations (9) and (10):
    ///
    /// ```text
    /// μ  = 2y / (2 + 2ax² − y²)
    /// x₃ = μx
    /// y₃ = μ(μ − y) − 1.
    /// ```
    ///
    pub fn double(&self, curve: &JacobiQuarticCurve<F>) -> Self {
        if self.is_identity() {
            return *self;
        }

        let x2 = <F as FieldOps>::square(&self.x);
        let y2 = <F as FieldOps>::square(&self.y);
        let two = <F as FieldOps>::double(&F::one());

        let denom = two + two * curve.a * x2 - y2;
        let denom_inv = denom.invert().into_option()
            .expect("Jacobi quartic doubling denominator vanished; result leaves affine chart or input is exceptional");

        let mu = two * self.y * denom_inv;
        let x3 = mu * self.x;
        let y3 = mu * (mu - self.y) - F::one();

        Self::new(x3, y3)
    }

    /// Affine addition using equations (1) and (2):
    ///
    /// ```text
    /// x₃ = (x₁y₂ + y₁x₂)/(1 − d x₁²x₂²)
    /// y₃ = ((y₁y₂ + 2ax₁x₂)(1 + d x₁²x₂²) + 2d x₁x₂(x₁² + x₂²))
    ///      /(1 − d x₁²x₂²)².
    /// ```
    ///
    pub fn add(&self, other: &Self, curve: &JacobiQuarticCurve<F>) -> Self {
        if self.is_identity() {
            return *other;
        }
        if other.is_identity() {
            return *self;
        }
        if *self == *other {
            return self.double(curve);
        }
        if *other == self.negate(curve) {
            return Self::identity();
        }

        let x1 = self.x;
        let y1 = self.y;
        let x2 = other.x;
        let y2 = other.y;

        let x1_sq = <F as FieldOps>::square(&x1);
        let x2_sq = <F as FieldOps>::square(&x2);
        let x1x2 = x1 * x2;
        let y1y2 = y1 * y2;
        let dx4 = curve.d * x1_sq * x2_sq;
        let two = <F as FieldOps>::double(&F::one());

        let denom = F::one() - dx4;
        let denom_inv = denom.invert().into_option()
            .expect("Jacobi quartic addition denominator vanished; choose d nonsquare / odd-order subgroup or use a projective model");

        let x3 = (x1 * y2 + y1 * x2) * denom_inv;

        let numer_y = (y1y2 + two * curve.a * x1x2) * (F::one() + dx4)
            + two * curve.d * x1x2 * (x1_sq + x2_sq);
        let y3 = numer_y * <F as FieldOps>::square(&denom_inv);

        Self::new(x3, y3)
    }

    /// Constant-time double-and-add in the same style as the Edwards code.
    pub fn scalar_mul(&self, k: &[u64], curve: &JacobiQuarticCurve<F>) -> Self {
        let mut result = Self::identity();

        for &limb in k.iter().rev() {
            for bit in (0..64).rev() {
                let doubled = result.double(curve);
                let added = doubled.add(self, curve);
                let choice = Choice::from(((limb >> bit) & 1) as u8);
                result = Self::conditional_select(&doubled, &added, choice);
            }
        }

        result
    }
}

impl<F: FieldOps> PointOps for JacobiQuarticPoint<F> {
    type BaseField = F;
    type Curve = JacobiQuarticCurve<F>;

    fn identity(_curve: &Self::Curve) -> Self {
        JacobiQuarticPoint::<F>::identity()
    }

    fn is_identity(&self) -> bool {
        JacobiQuarticPoint::<F>::is_identity(self)
    }

    fn negate(&self, curve: &Self::Curve) -> Self {
        JacobiQuarticPoint::<F>::negate(self, curve)
    }

    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        JacobiQuarticPoint::<F>::scalar_mul(self, k, curve)
    }
}

impl<F: FieldOps> PointAdd for JacobiQuarticPoint<F> {
    fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
        JacobiQuarticPoint::<F>::add(self, other, curve)
    }
}
