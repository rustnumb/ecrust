//! Point representation and group law for Edwards curves.
//!
//! A single `EdwardsPoint<F>` type handles both odd and even characteristic
//! via runtime dispatch on `F::characteristic()`.
//!
//! # Odd characteristic (`x² + y² = 1 + d·x²·y²`)
//!
//! - Identity: `(0, 1)`
//! - Negation: `-(x, y) = (-x, y)`
//! - Addition: `x₃ = (x₁y₂+y₁x₂)/(1+d·x₁x₂y₁y₂)`,
//!             `y₃ = (y₁y₂-x₁x₂)/(1-d·x₁x₂y₁y₂)`
//!
//! # Characteristic 2 (`d₁(x+y)+d₂(x²+y²) = xy+xy(x+y)+x²y²`)
//!
//! - Identity: `(0, 0)`
//! - Negation: `-(x, y) = (y, x)`
//! - Addition: Bernstein–Lange–Rezaeian Farashahi §3 formulas
//!   (strongly unified — works for doubling too)

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_edwards::EdwardsCurve;
use crate::point_ops::PointOps;
use fp::field_ops::FieldOps;

/// An affine point on an Edwards curve, for any characteristic.
#[derive(Debug, Clone, Copy)]
pub struct EdwardsPoint<F: FieldOps> {
    pub x: F,
    pub y: F,
}

// ---------------------------------------------------------------------------
// PartialEq / Eq
// ---------------------------------------------------------------------------

impl<F: FieldOps> PartialEq for EdwardsPoint<F> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<F: FieldOps> Eq for EdwardsPoint<F> {}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<F: FieldOps> EdwardsPoint<F> {
    /// Construct an affine Edwards point.  No on-curve check.
    pub fn new(x: F, y: F) -> Self {
        Self { x, y }
    }

    /// The group identity.
    /// - Odd char: `(0, 1)`
    /// - Char 2:   `(0, 0)`
    pub fn identity() -> Self {
        if F::characteristic()[0] != 2 {
            Self { x: F::zero(), y: F::one() }
        } else {
            Self { x: F::zero(), y: F::zero() }
        }
    }

    /// Check whether this is the identity element.
    pub fn is_identity(&self) -> bool {
        if F::characteristic()[0] != 2 {
            self.x == F::zero() && self.y == F::one()
        } else {
            self.x == F::zero() && self.y == F::zero()
        }
    }
}

// ---------------------------------------------------------------------------
// Constant-time
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for EdwardsPoint<F>
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

impl<F> ConstantTimeEq for EdwardsPoint<F>
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

// ---------------------------------------------------------------------------
// Group operations
// ---------------------------------------------------------------------------

impl<F: FieldOps> EdwardsPoint<F> {
    /// Negate a point.
    /// - Odd char: `-(x, y) = (-x, y)`
    /// - Char 2:   `-(x, y) = (y, x)`
    pub fn negate(&self, _curve: &EdwardsCurve<F>) -> Self {
        if F::characteristic()[0] != 2 {
            Self::new(-self.x, self.y)
        } else {
            Self::new(self.y, self.x)
        }
    }

    /// Add two points on the Edwards curve.
    pub fn add(&self, other: &Self, curve: &EdwardsCurve<F>) -> Self {
        if F::characteristic()[0] != 2 {
            self.add_odd(other, curve)
        } else {
            self.add_binary(other, curve)
        }
    }

    /// Double a point.  Both addition laws are strongly unified, so this
    /// just delegates to `add`.
    pub fn double(&self, curve: &EdwardsCurve<F>) -> Self {
        self.add(self, curve)
    }

    /// Scalar multiplication `[k]P` (constant-time double-and-add).
    pub fn scalar_mul(&self, k: &[u64], curve: &EdwardsCurve<F>) -> Self {
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

    // -----------------------------------------------------------------------
    // Odd-characteristic addition
    //
    //   x₃ = (x₁y₂ + y₁x₂) / (1 + d·x₁x₂y₁y₂)
    //   y₃ = (y₁y₂ − x₁x₂) / (1 − d·x₁x₂y₁y₂)
    // -----------------------------------------------------------------------

    fn add_odd(&self, other: &Self, curve: &EdwardsCurve<F>) -> Self {
        let d = curve.d2;

        let x1y2 = self.x * other.y;
        let y1x2 = self.y * other.x;
        let x1x2 = self.x * other.x;
        let y1y2 = self.y * other.y;

        let dxy = d * x1x2 * y1y2;

        let one = F::one();
        let x_num = x1y2 + y1x2;
        let x_den = one + dxy;
        let y_num = y1y2 - x1x2;
        let y_den = one - dxy;

        let x_den_inv = x_den.invert().into_option()
            .expect("Edwards addition: x-denominator must be invertible");
        let y_den_inv = y_den.invert().into_option()
            .expect("Edwards addition: y-denominator must be invertible");

        Self::new(x_num * x_den_inv, y_num * y_den_inv)
    }

    // -----------------------------------------------------------------------
    // Characteristic-2 addition  (Bernstein–Lange–Rezaeian Farashahi §3)
    //
    //   w₁ = x₁+y₁,  w₂ = x₂+y₂
    //   A  = x₁²+x₁,  B  = y₁²+y₁
    //
    //   x₃ = (d₁(x₁+x₂) + d₂·w₁·w₂ + A·(x₂(y₁+y₂+1)+y₁y₂))
    //         / (d₁ + A·w₂)
    //   y₃ = (d₁(y₁+y₂) + d₂·w₁·w₂ + B·(y₂(x₁+x₂+1)+x₁x₂))
    //         / (d₁ + B·w₂)
    // -----------------------------------------------------------------------

    fn add_binary(&self, other: &Self, curve: &EdwardsCurve<F>) -> Self {
        let x1 = self.x;
        let y1 = self.y;
        let x2 = other.x;
        let y2 = other.y;
        let d1 = curve.d1;
        let d2 = curve.d2;

        let w1 = x1 + y1;
        let w2 = x2 + y2;

        let a = <F as FieldOps>::square(&x1) + x1;
        let b = <F as FieldOps>::square(&y1) + y1;

        let d2_w1w2 = d2 * w1 * w2;
        let x1x2 = x1 * x2;
        let y1y2 = y1 * y2;

        let x_num = d1 * (x1 + x2) + d2_w1w2
            + a * (x2 * (y1 + y2 + F::one()) + y1y2);
        let x_den = d1 + a * w2;

        let y_num = d1 * (y1 + y2) + d2_w1w2
            + b * (y2 * (x1 + x2 + F::one()) + x1x2);
        let y_den = d1 + b * w2;

        let x_den_inv = x_den.invert().into_option()
            .expect("binary Edwards addition: x-denom must be invertible");
        let y_den_inv = y_den.invert().into_option()
            .expect("binary Edwards addition: y-denom must be invertible");

        Self::new(x_num * x_den_inv, y_num * y_den_inv)
    }
}

// ---------------------------------------------------------------------------
// w-coordinate helpers for char-2 Montgomery ladder
// ---------------------------------------------------------------------------

impl<F: FieldOps> EdwardsPoint<F> {
    /// Differential addition on the `w`-line (`w = x + y`, char 2 only).
    ///
    /// Given `w₁ = w(Q−P)`, `w₂ = w(P)`, `w₃ = w(Q)`, compute `w₅ = w(P+Q)`.
    pub fn w_diff_add(w1: &F, w2: &F, w3: &F, curve: &EdwardsCurve<F>) -> F {
        let r = *w2 * *w3;
        let s = <F as FieldOps>::square(&r);
        let one = F::one();
        let t = r * (one + *w2 + *w3) + s;

        let d2_over_d1 = curve.d2
            * curve.d1.invert().into_option().expect("d1 invertible");
        let coeff = d2_over_d1 + one;

        let den = curve.d1 + t + coeff * s;
        let den_inv = den.invert().into_option()
            .expect("w-diff-add denominator invertible");

        t * den_inv + *w1
    }

    /// `w`-coordinate doubling (`w = x + y`, char 2 only).
    ///
    /// Given `w₂ = w(P)`, compute `w₄ = w(2P)`.
    pub fn w_double(w2: &F, curve: &EdwardsCurve<F>) -> F {
        let a = <F as FieldOps>::square(w2);
        let j = <F as FieldOps>::square(&a);

        let d2_over_d1 = curve.d2
            * curve.d1.invert().into_option().expect("d1 invertible");

        let num = a + j;
        let den = curve.d1 + a + d2_over_d1 * j;
        let den_inv = den.invert().into_option()
            .expect("w-double denominator invertible");

        num * den_inv
    }
}

// ---------------------------------------------------------------------------
// PointOps bridge
// ---------------------------------------------------------------------------

impl<F: FieldOps> PointOps for EdwardsPoint<F> {
    type BaseField = F;
    type Curve = EdwardsCurve<F>;

    fn identity(_curve: &Self::Curve) -> Self {
        EdwardsPoint::<F>::identity()
    }

    fn is_identity(&self) -> bool {
        EdwardsPoint::<F>::is_identity(self)
    }

    fn negate(&self, curve: &Self::Curve) -> Self {
        EdwardsPoint::<F>::negate(self, curve)
    }

    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        EdwardsPoint::<F>::scalar_mul(self, k, curve)
    }
}
