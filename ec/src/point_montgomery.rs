//! X-only arithmetic on a Montgomery curve via the Kummer line.
//!
//!
//! # Kummer-line representation
//!
//! A Montgomery point is represented only by its x-coordinate, in projective
//! form `(X : Z)`, corresponding to the affine value `x = X / Z` when `Z ≠ 0`.
//!
//! This is not a full point on the elliptic curve: it is a point on the
//! quotient `E / {±1}`. In other words, `P` and `-P` are identified.
//!
//! # Why projective x/z coordinates?
//!
//! The Montgomery ladder is naturally expressed using x-only formulas in
//! projective coordinates. This avoids field inversions during the ladder and
//! yields a uniform constant-time scalar multiplication pattern.
//!
//! # Important limitation
//!
//! Since the sign of `y` is forgotten, x-only arithmetic does **not** provide a
//! full group law on points. In particular, ordinary addition `P + Q` is not
//! available from x-coordinates alone.
//!
//! What *is* available is:
//!
//! - doubling,
//! - differential addition: given `x(P)`, `x(Q)`, and `x(P - Q)`, compute
//!   `x(P + Q)`,
//! - scalar multiplication via the Montgomery ladder.

use subtle::{Choice, CtOption, ConditionallySelectable, ConstantTimeEq};

use crate::curve_montgomery::MontgomeryCurve;
use crate::point_ops::PointOps;
use fp::field_ops::FieldOps;

/// A point on the Kummer line of a Montgomery curve, represented by `(X : Z)`.
///
/// The affine x-coordinate is `x = X / Z` when `Z ≠ 0`.
/// A conventional choice is:
/// - `(1 : 0)` for the identity image,
/// - `(X : Z)` with `Z ≠ 0` for finite x-coordinates.
#[derive(Debug, Clone, Copy)]
pub struct KummerPoint<F: FieldOps> {
    pub x: F,
    pub z: F,
}

// ---------------------------------------------------------------------------
// Manual trait impls
// ---------------------------------------------------------------------------

impl<F: FieldOps> PartialEq for KummerPoint<F>
where
    F: FieldOps + ConstantTimeEq,
{
    /// Equality of projective x-line points.
    ///
    /// A standard criterion is cross-multiplication:
    /// ```text
    /// X1 Z2 = X2 Z1.
    /// ```
    fn eq(&self, other: &Self) -> bool {
        self.x * other.z == other.x * self.z
    }
}

impl<F: FieldOps> Eq for KummerPoint<F>
where
    F: FieldOps + ConstantTimeEq,
{ }

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<F: FieldOps> KummerPoint<F> {
    /// Construct a projective x-line point `(X : Z)` without validation.
    pub fn new(x: F, z: F) -> Self {
        Self{ x, z }
    }

    /// Construct the finite x-line point corresponding to the affine
    /// x-coordinate `x`, i.e. `(x : 1)`.
    pub fn from_x(x: F) -> Self {
        Self{ x, z: F::one() }
    }

    /// The image of the identity point on the Kummer line.
    pub fn identity() -> Self {
        Self{ x: F::one(), z: F::zero()}
    }

    /// Return `true` if this point is the image of identity.
    pub fn is_identity(&self) -> bool {
        bool::from(self.z.is_zero())
    }

    /// Attempt to recover the affine x-coordinate (succeeds when Z != 0).
    pub fn to_x(&self) -> CtOption<F> {
        let z_is_non_zero = !self.z.is_zero();
        CtOption::new(self.x, z_is_non_zero)
    }
}


// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for KummerPoint<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self{
            x: F::conditional_select(& a.x, &b.x, choice),
            z: F::conditional_select(& a.z, &b.z, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.x, &other.x, choice);
        F::conditional_assign(&mut self.z, &other.z, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.x, &mut b.x, choice);
        F::conditional_swap(&mut a.z, &mut b.z, choice);
    }
}

impl<F> ConstantTimeEq for KummerPoint<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    /// Constant-time projective equality test on the Kummer line.
    /// ```text
    /// X1 Z2  ?=  X2 Z1
    /// ```
    fn ct_eq(&self, other: &Self) -> Choice {
        let x1z2 = self.x * other.z;
        let x2z1 = other.x * self.z;
        x1z2.ct_eq(& x2z1)
    }

    fn ct_ne(&self, other: &Self) -> Choice { !self.ct_eq(other) }
}


// ---------------------------------------------------------------------------
// X-only arithmetic
// ---------------------------------------------------------------------------

impl<F: FieldOps> KummerPoint<F> {
    /// Point doubling on the Kummer line.
    ///
    /// Given `x(P)` in projective form `(X:Z)`, compute `x([2]P)`.
    pub fn xdouble(&self, curve: &MontgomeryCurve<F>) -> Self {
        if F::characteristic()[0] != 2 {
            let q = <F as FieldOps>::square(&(self.x + self.z));
            let r = <F as FieldOps>::square(&(self.x - self.z));
            let s = q - r;
            let a24 = curve.a24();

            let new_z = s * (r + a24 * s);

            Self{ x: q * r, z: new_z }
        }
        else {
            let temp_x = self.x + curve.b * self.z;
            let new_x = <F as FieldOps>::square(&<F as FieldOps>::square(&temp_x));
            let xz = self.x * self.z;
            let new_z = <F as FieldOps>::square(&xz);
            Self{ x: new_x, z: new_z }
        }

    }

    /// Differential addition. Given (in projective form `(X:Z)`):
    ///
    /// - `self = x(P)`,
    /// - `other = x(Q)`,
    /// - `diff = x(P - Q)`,
    ///
    /// compute `x(P + Q)`.
    pub fn xadd(&self, other: &Self, diff: &Self) -> Self {
        if F::characteristic()[0] != 2 {
            let u = (self.x - self.z) * (other.x + other.z);
            let v = (self.x + self.z) * (other.x - other.z);

            let new_x = diff.z * <F as FieldOps>::square(&(u + v));
            let new_z = diff.x * <F as FieldOps>::square(&(u - v));

            Self { x: new_x, z: new_z }
        }
        else {
            let x1z2 = self.x * other.z;
            let x2z1 = other.x * self.z;
            let sum_squared = <F as FieldOps>::square(&(x1z2 + x2z1));

            let new_x = diff.x * sum_squared + diff.z * x1z2 * x2z1;
            let new_z = diff.z * sum_squared;

            Self { x: new_x, z: new_z }
        }
    }

    /// Montgomery ladder for scalar multiplication.
    ///
    /// Given an x-line point `x(P)` and a scalar `k`, compute `x([k]P)`.
    /// The scalar `k` is given as a slice of `u64` limbs in **little-endian**
    /// order (same convention as `FieldOps::pow`).
    pub fn scalar_mul(&self, k: &[u64], curve: &MontgomeryCurve<F>) -> Self {
        if self.is_identity() || k.is_empty() {
            return Self::identity();
        }

        let mut result = Self::identity();
        let diff = self.clone();

        for &limb in k.iter().rev() {
            for bit in (0..64).rev() {
                let doubled = result.xdouble(curve);
                let added = doubled.xadd(self, &diff);
                let choice = Choice::from(((limb >> bit) & 1) as u8);
                result = Self::conditional_select(&doubled, &added, choice);
            }
        }

        result
    }
}

// ---------------------------------------------------------------------------
// PointOps bridge
// ---------------------------------------------------------------------------

impl<F> PointOps for KummerPoint<F>
where
    F: FieldOps,
{
    type BaseField = F;
    type Curve = MontgomeryCurve<F>;

    /// Return the identity image on the Kummer line.
    fn identity(_curve: &Self::Curve) -> Self {
        KummerPoint::<F>::identity()
    }

    /// Return `true` if this is the identity image.
    fn is_identity(&self) -> bool {
        KummerPoint::<F>::is_identity(self)
    }

    /// Negation is trivial on the Kummer line because `P` and `-P` have the
    /// same image.
    fn negate(&self, _curve: &Self::Curve) -> Self {
        *self
    }

    /// Scalar multiplication is naturally implemented by the Montgomery ladder.
    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        KummerPoint::<F>::scalar_mul(self, k, curve)
    }
}