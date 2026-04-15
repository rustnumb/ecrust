//! Affine points on a Jacobi-intersection curve.
//!
//! The EFD gives the affine group law for
//!
//! $$
//! s^2 + c^2 = 1 \quad
//! a s^2 + d^2 = 1.
//! $$
//!
//! We use those formulas directly, with neutral element `(0, 1, 1)` and
//! negation `-(s, c, d) = (-s, c, d)`.

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_jacobi_intersection::JacobiIntersectionCurve;
use crate::point_ops::{PointAdd, PointOps};
use fp::field_ops::FieldOps;

/// An affine point $(s, c, d)$ on a Jacobi intersection.
#[derive(Debug, Clone, Copy)]
pub struct JacobiIntersectionPoint<F: FieldOps> {
    /// The coordinate `s` of the point
    pub s: F,
    /// The coordinate `c` of the point
    pub c: F,
    /// The coordinate `d` of the point
    pub d: F,
}

impl<F: FieldOps> PartialEq for JacobiIntersectionPoint<F> {
    fn eq(&self, other: &Self) -> bool {
        self.s == other.s && self.c == other.c && self.d == other.d
    }
}

impl<F: FieldOps> Eq for JacobiIntersectionPoint<F> {}

impl<F: FieldOps> JacobiIntersectionPoint<F> {
    /// Creates a new point on a Jacobi intersection
    ///
    /// The point should satisfy
    /// $$
    /// s^2 + c^2 = 1 \quad
    /// a s^2 + d^2 = 1.
    /// $$
    ///
    /// # Arguments
    ///
    /// * `s` - Element of `F` (type: F)
    /// * `c` - Element of `F` (type: F)
    /// * `d` - Element of `F` (type: F)
    ///
    /// # Returns
    ///
    /// Returns the point on the Jacobi intersection (type: Self)
    pub fn new(s: F, c: F, d: F) -> Self {
        Self { s, c, d }
    }

    /// The neutral element $(0, 1, 1)$.
    pub fn identity() -> Self {
        Self {
            s: F::zero(),
            c: F::one(),
            d: F::one(),
        }
    }

    /// Checks if the point is the indentity point
    ///
    /// # Arguments
    ///
    /// * `&self` - A point on the curve (type: self)
    ///
    /// # Returns
    ///
    /// True if and only if it is $(0,1,1)$ (type: bool)
    pub fn is_identity(&self) -> bool {
        self.s == F::zero() && self.c == F::one() && self.d == F::one()
    }
}

impl<F> ConditionallySelectable for JacobiIntersectionPoint<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            s: F::conditional_select(&a.s, &b.s, choice),
            c: F::conditional_select(&a.c, &b.c, choice),
            d: F::conditional_select(&a.d, &b.d, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.s.conditional_assign(&other.s, choice);
        self.c.conditional_assign(&other.c, choice);
        self.d.conditional_assign(&other.d, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.s, &mut b.s, choice);
        F::conditional_swap(&mut a.c, &mut b.c, choice);
        F::conditional_swap(&mut a.d, &mut b.d, choice);
    }
}

impl<F> ConstantTimeEq for JacobiIntersectionPoint<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.s.ct_eq(&other.s) & self.c.ct_eq(&other.c) & self.d.ct_eq(&other.d)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}

impl<F: FieldOps> JacobiIntersectionPoint<F> {
    /// Negation: `-(s, c, d) = (-s, c, d)`.
    pub fn negate(&self, _curve: &JacobiIntersectionCurve<F>) -> Self {
        Self::new(-self.s, self.c, self.d)
    }

    /// Affine addition formulas from the EFD:
    ///
    /// ```text
    /// s₃ = (c₂ s₁ d₂ + d₁ s₂ c₁)/(c₂² + (d₁ s₂)²)
    /// c₃ = (c₂ c₁ - d₁ s₂ s₁ d₂)/(c₂² + (d₁ s₂)²)
    /// d₃ = (d₁ d₂ - a s₁ c₁ s₂ c₂)/(c₂² + (d₁ s₂)²).
    /// ```
    ///
    pub fn add(&self, other: &Self, curve: &JacobiIntersectionCurve<F>) -> Self {
        let d1s2 = self.d * other.s;
        let denom = <F as FieldOps>::square(&other.c) + <F as FieldOps>::square(&d1s2);
        let denom_inv = denom
            .invert()
            .into_option()
            .expect("Jacobi-intersection addition denominator vanished");

        let s3 = (other.c * self.s * other.d + self.d * other.s * self.c) * denom_inv;
        let c3 = (other.c * self.c - self.d * other.s * self.s * other.d) * denom_inv;
        let d3 = (self.d * other.d - curve.a * self.s * self.c * other.s * other.c) * denom_inv;

        Self::new(s3, c3, d3)
    }

    /// Affine doubling formulas from the EFD:
    ///
    /// ```text
    /// s₃ = 2c s d / (c² + (d s)²)
    /// c₃ = (c² - d² s²) / (c² + (d s)²)
    /// d₃ = (d² - a s² c²) / (c² + (d s)²).
    /// ```
    ///
    pub fn double(&self, curve: &JacobiIntersectionCurve<F>) -> Self {
        let ds = self.d * self.s;
        let denom = <F as FieldOps>::square(&self.c) + <F as FieldOps>::square(&ds);
        let denom_inv = denom
            .invert()
            .into_option()
            .expect("Jacobi-intersection doubling denominator vanished");

        let two = <F as FieldOps>::double(&F::one());
        let s_sq = <F as FieldOps>::square(&self.s);
        let c_sq = <F as FieldOps>::square(&self.c);
        let d_sq = <F as FieldOps>::square(&self.d);

        let s3 = (two * self.c * self.s * self.d) * denom_inv;
        let c3 = (c_sq - d_sq * s_sq) * denom_inv;
        let d3 = (d_sq - curve.a * s_sq * c_sq) * denom_inv;

        Self::new(s3, c3, d3)
    }

    /// Scalar multiplication on the curve
    ///
    /// Multiplies a point by an integers
    ///
    /// # Arguments
    ///
    /// * `&self` - Point on the curve (type: `Self`)
    /// * `k` - An integer represented as a vec of `u64` (type: `&[u64]`)
    /// * `curve` - The curve on which `self` lies (type: `&JacobiIntersectionCurve<F>`)
    ///
    /// # Returns
    ///
    /// The point `k * self` (type: Self)
    pub fn scalar_mul(&self, k: &[u64], curve: &JacobiIntersectionCurve<F>) -> Self {
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

impl<F: FieldOps> PointOps for JacobiIntersectionPoint<F> {
    type BaseField = F;
    type Curve = JacobiIntersectionCurve<F>;

    fn identity(_curve: &Self::Curve) -> Self {
        JacobiIntersectionPoint::<F>::identity()
    }

    fn is_identity(&self) -> bool {
        JacobiIntersectionPoint::<F>::is_identity(self)
    }

    fn negate(&self, curve: &Self::Curve) -> Self {
        JacobiIntersectionPoint::<F>::negate(self, curve)
    }

    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        JacobiIntersectionPoint::<F>::scalar_mul(self, k, curve)
    }
}

impl<F: FieldOps> PointAdd for JacobiIntersectionPoint<F> {
    fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
        JacobiIntersectionPoint::<F>::add(self, other, curve)
    }
}
