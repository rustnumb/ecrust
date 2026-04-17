//! Projective points on a generalized Hessian curve.
//!
//! We represent points on
//!
//! $$X^3 + Y^3 + cZ^3 = dXYZ$$
//!
//! by projective triples `(X:Y:Z)`.
//!
//! This choice is important for Hessian curves because the neutral element is a
//! point at infinity:
//!
//! $$O = (1 : -1 : 0).$$
//!
//! The formulas implemented here follow:
//!
//! - Farashahi--Joye, §2--§4 for the generalized Hessian model,
//! - the EFD projective Hessian formulas for the ordinary Hessian case.
//!
//! In particular:
//!
//! - negation is `-(X:Y:Z) = (Y:X:Z)`,
//! - doubling uses the projective formulas from equation (6) in the paper,
//! - addition uses the unified formulas (9), with formulas (10) as a fallback
//!   for the exceptional cases described in §4.

use core::fmt;

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_hessian::HessianCurve;
use crate::point_ops::{PointAdd, PointOps};
use fp::field_ops::FieldOps;

/// A projective point `(X:Y:Z)` on a generalized Hessian curve.
#[derive(Debug, Clone, Copy)]
pub struct HessianPoint<F: FieldOps> {
    /// Projective `X` coordinate.
    pub x: F,
    /// Projective `Y` coordinate.
    pub y: F,
    /// Projective `Z` coordinate.
    pub z: F,
}

impl<F> fmt::Display for HessianPoint<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_identity() {
            if f.alternate() {
                write!(f, "HessianPoint {{ O = (1:-1:0) }}")
            } else {
                write!(f, "O")
            }
        } else if self.is_zero_projective() {
            if f.alternate() {
                write!(f, "HessianPoint {{ invalid = (0:0:0) }}")
            } else {
                write!(f, "(0:0:0)")
            }
        } else if let Some((x_aff, y_aff)) = self.to_affine() {
            if f.alternate() {
                write!(
                    f,
                    "HessianPoint {{ X:Y:Z = ({}:{}:{}), x = {}, y = {} }}",
                    self.x, self.y, self.z, x_aff, y_aff
                )
            } else {
                write!(f, "({}, {})", x_aff, y_aff)
            }
        } else if f.alternate() {
            write!(f, "HessianPoint {{ X:Y:Z = ({}:{}:{}) }}", self.x, self.y, self.z)
        } else {
            write!(f, "({}:{}:{})", self.x, self.y, self.z)
        }
    }
}

impl<F: FieldOps> PartialEq for HessianPoint<F> {
    /// Equality of projective points.
    fn eq(&self, other: &Self) -> bool {
        let self_zero = self.is_zero_projective();
        let other_zero = other.is_zero_projective();
        if self_zero || other_zero {
            return self_zero && other_zero;
        }

        self.x * other.y == other.x * self.y
            && self.x * other.z == other.x * self.z
            && self.y * other.z == other.y * self.z
    }
}

impl<F: FieldOps> Eq for HessianPoint<F> {}

impl<F: FieldOps> HessianPoint<F> {
    /// Construct a projective Hessian point without validation.
    pub fn new(x: F, y: F, z: F) -> Self {
        Self { x, y, z }
    }

    /// Construct the finite affine point `(x, y)`, represented as `(x:y:1)`.
    pub fn from_affine(x: F, y: F) -> Self {
        Self { x, y, z: F::one() }
    }

    /// Return the neutral element `(1:-1:0)`.
    pub fn identity() -> Self {
        Self {
            x: F::one(),
            y: -F::one(),
            z: F::zero(),
        }
    }

    /// Return `true` when the point is the neutral element.
    pub fn is_identity(&self) -> bool {
        if !bool::from(self.z.is_zero()) {
            return false;
        }

        if bool::from(self.x.is_zero()) && bool::from(self.y.is_zero()) {
            return false;
        }

        self.x + self.y == F::zero()
    }

    /// Return `true` if the point lies on the line at infinity.
    pub fn is_at_infinity(&self) -> bool {
        bool::from(self.z.is_zero()) && !self.is_zero_projective()
    }

    /// Return `true` if this is the invalid projective triple `(0:0:0)`.
    pub fn is_zero_projective(&self) -> bool {
        bool::from(self.x.is_zero())
            && bool::from(self.y.is_zero())
            && bool::from(self.z.is_zero())
    }

    /// Convert a finite projective point to affine coordinates.
    pub fn to_affine(&self) -> Option<(F, F)> {
        self.z
            .invert()
            .into_option()
            .map(|zinv| (self.x * zinv, self.y * zinv))
    }
}

impl<F> ConditionallySelectable for HessianPoint<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            x: F::conditional_select(&a.x, &b.x, choice),
            y: F::conditional_select(&a.y, &b.y, choice),
            z: F::conditional_select(&a.z, &b.z, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.x.conditional_assign(&other.x, choice);
        self.y.conditional_assign(&other.y, choice);
        self.z.conditional_assign(&other.z, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.x, &mut b.x, choice);
        F::conditional_swap(&mut a.y, &mut b.y, choice);
        F::conditional_swap(&mut a.z, &mut b.z, choice);
    }
}

impl<F> ConstantTimeEq for HessianPoint<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        let self_zero = self.is_zero_projective();
        let other_zero = other.is_zero_projective();
        if self_zero || other_zero {
            return Choice::from((self_zero && other_zero) as u8);
        }

        (self.x * other.y).ct_eq(&(other.x * self.y))
            & (self.x * other.z).ct_eq(&(other.x * self.z))
            & (self.y * other.z).ct_eq(&(other.y * self.z))
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}

impl<F: FieldOps> HessianPoint<F> {
    /// Negation on a Hessian curve:
    ///
    /// $$-(X:Y:Z) = (Y:X:Z).$$
    pub fn negate(&self, _curve: &HessianCurve<F>) -> Self {
        Self::new(self.y, self.x, self.z)
    }

    /// Projective point doubling.
    ///
    /// This implements Farashahi--Joye equation (6):
    ///
    /// ```text
    /// X3 = Y1 (c Z1^3 - X1^3)
    /// Y3 = X1 (Y1^3 - c Z1^3)
    /// Z3 = Z1 (X1^3 - Y1^3)
    /// ```
    pub fn double(&self, curve: &HessianCurve<F>) -> Self {
        if self.is_identity() {
            return *self;
        }

        let x2 = <F as FieldOps>::square(&self.x);
        let y2 = <F as FieldOps>::square(&self.y);
        let z2 = <F as FieldOps>::square(&self.z);

        let x3 = self.x * x2;
        let y3 = self.y * y2;
        let z3 = self.z * z2;
        let cz3 = curve.c * z3;

        Self::new(
            self.y * (cz3 - x3),
            self.x * (y3 - cz3),
            self.z * (x3 - y3),
        )
    }

    /// Unified projective addition formula (9) from Farashahi--Joye §3.
    fn add_formula_9(&self, other: &Self, curve: &HessianCurve<F>) -> Self {
        let x1_sq = <F as FieldOps>::square(&self.x);
        let y1_sq = <F as FieldOps>::square(&self.y);
        let z1_sq = <F as FieldOps>::square(&self.z);
        let x2_sq = <F as FieldOps>::square(&other.x);
        let y2_sq = <F as FieldOps>::square(&other.y);
        let z2_sq = <F as FieldOps>::square(&other.z);

        let x3 = curve.c * other.y * other.z * z1_sq - self.x * self.y * x2_sq;
        let y3 = other.x * other.y * y1_sq - curve.c * self.x * self.z * z2_sq;
        let z3 = other.x * other.z * x1_sq - self.y * self.z * y2_sq;

        Self::new(x3, y3, z3)
    }

    /// Unified projective addition formula (10) from Farashahi--Joye §3.
    fn add_formula_10(&self, other: &Self, curve: &HessianCurve<F>) -> Self {
        let x1_sq = <F as FieldOps>::square(&self.x);
        let y1_sq = <F as FieldOps>::square(&self.y);
        let z1_sq = <F as FieldOps>::square(&self.z);
        let x2_sq = <F as FieldOps>::square(&other.x);
        let y2_sq = <F as FieldOps>::square(&other.y);
        let z2_sq = <F as FieldOps>::square(&other.z);

        let x3 = curve.c * self.y * self.z * z2_sq - other.x * other.y * x1_sq;
        let y3 = self.x * self.y * y2_sq - curve.c * other.x * other.z * z1_sq;
        let z3 = self.x * self.z * x2_sq - other.y * other.z * y1_sq;

        Self::new(x3, y3, z3)
    }

    /// Add two projective Hessian points.
    ///
    /// We first evaluate the unified formulas (9). If they produce the invalid
    /// triple `(0:0:0)`, we fall back to formulas (10), which cover the
    /// complementary exceptional set described in Farashahi--Joye §4.
    pub fn add(&self, other: &Self, curve: &HessianCurve<F>) -> Self {
        if self.is_identity() {
            return *other;
        }
        if other.is_identity() {
            return *self;
        }

        let r = self.add_formula_9(other, curve);
        if !r.is_zero_projective() {
            return r;
        }

        let s = self.add_formula_10(other, curve);
        if !s.is_zero_projective() {
            return s;
        }

        if *other == self.negate(curve) {
            return Self::identity();
        }

        panic!("Hessian addition failed for valid-looking inputs; both unified formula branches vanished");
    }

    /// Variable-time double-and-add scalar multiplication.
    pub fn scalar_mul(&self, k: &[u64], curve: &HessianCurve<F>) -> Self {
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

impl<F: FieldOps> PointOps for HessianPoint<F> {
    type BaseField = F;
    type Curve = HessianCurve<F>;

    fn identity(_curve: &Self::Curve) -> Self {
        HessianPoint::<F>::identity()
    }

    fn is_identity(&self) -> bool {
        HessianPoint::<F>::is_identity(self)
    }

    fn negate(&self, curve: &Self::Curve) -> Self {
        HessianPoint::<F>::negate(self, curve)
    }

    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        HessianPoint::<F>::scalar_mul(self, k, curve)
    }
}

impl<F: FieldOps> PointAdd for HessianPoint<F> {
    fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
        HessianPoint::<F>::add(self, other, curve)
    }
}
