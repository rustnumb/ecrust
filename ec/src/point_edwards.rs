//! Point representation and group law for Edwards curves.
//!
//! A single `EdwardsPoint<F>` type handles both odd and even characteristic
//! via runtime dispatch on $F::characteristic()$.
//!
//! # Odd characteristic
//!
//! The curve is given by
//!
//! $$
//! x^2 + y^2 = 1 + d x^2 y^2
//! $$
//!
//! - Identity: $(0, 1)$
//! - Negation: $-(x, y) = (-x, y)$
//! - Addition:
//!
//! $$
//! x_3 = \frac{x_1 y_2 + y_1 x_2}{1 + d x_1 x_2 y_1 y_2},
//! \quad
//! y_3 = \frac{y_1 y_2 - x_1 x_2}{1 - d x_1 x_2 y_1 y_2}
//! $$
//!
//! # Characteristic $2$
//!
//! The curve is given by
//!
//! $$
//! d_1(x + y) + d_2(x^2 + y^2) = xy + xy(x + y) + x^2 y^2
//! $$
//!
//! - Identity: $(0, 0)$
//! - Negation: $-(x, y) = (y, x)$
//! - Addition: Bernstein–Lange–Rezaeian Farashahi §3 formulas
//!   (strongly unified — works for doubling too)

use core::fmt;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_edwards::EdwardsCurve;
use crate::point_ops::PointOps;
use fp::field_ops::FieldOps;
use fp::ref_field_impl;

/// An affine point on an Edwards curve, for any characteristic.
#[derive(Debug, Clone, Copy)]
pub struct EdwardsPoint<F: FieldOps> {
    /// The x coordinate of a point.
    pub x: F,
    /// The y coordinate of a point.
    pub y: F,
}

impl<F> fmt::Display for EdwardsPoint<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_identity() {
            if f.alternate() {
                if F::characteristic()[0] != 2 {
                    write!(f, "EdwardsPoint {{ O = (0,1) }}")
                } else {
                    write!(f, "EdwardsPoint {{ O = (0,0) }}")
                }
            } else {
                write!(f, "O")
            }
        } else if f.alternate() {
            write!(f, "EdwardsPoint {{ x = {}, y = {} }}", self.x, self.y)
        } else {
            write!(f, "({}, {})", self.x, self.y)
        }
    }
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
    /// Construct an affine Edwards point. No on-curve check.
    pub fn new(x: F, y: F) -> Self {
        Self { x, y }
    }

    /// The group identity.
    /// - Odd char: `(0, 1)`
    /// - Char 2:   `(0, 0)`
    pub fn identity() -> Self {
        if F::characteristic()[0] != 2 {
            Self {
                x: F::zero(),
                y: F::one(),
            }
        } else {
            Self {
                x: F::zero(),
                y: F::zero(),
            }
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

ref_field_impl! {
    impl<F> EdwardsPoint<F> {
        /// Negate a point.
        /// - Odd char: `-(x, y) = (-x, y)`
        /// - Char 2:   `-(x, y) = (y, x)`
        pub fn negate(&self, _curve: &EdwardsCurve<F>) -> Self {
            if F::characteristic()[0] != 2 {
                Self::new(-&self.x, self.y.clone())
            } else {
                Self::new(self.y.clone(), self.x.clone())
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

        /// Double a point. Both addition laws are strongly unified, so this
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
            let x1y2 = &self.x * &other.y;
            let y1x2 = &self.y * &other.x;
            let x1x2 = &self.x * &other.x;
            let y1y2 = &self.y * &other.y;

            let x1x2y1y2 = &x1x2 * &y1y2;
            let dxy = &curve.d2 * &x1x2y1y2;

            let one = F::one();
            let x_num = &x1y2 + &y1x2;
            let x_den = &one + &dxy;
            let y_num = &y1y2 - &x1x2;
            let y_den = &one - &dxy;

            let x_den_inv = <F as FieldOps>::invert(&x_den)
                .into_option()
                .expect("Edwards addition: x-denominator must be invertible");
            let y_den_inv = <F as FieldOps>::invert(&y_den)
                .into_option()
                .expect("Edwards addition: y-denominator must be invertible");

            let x3 = &x_num * &x_den_inv;
            let y3 = &y_num * &y_den_inv;

            Self::new(x3, y3)
        }

        // -----------------------------------------------------------------------
        // Characteristic-2 addition
        //
        //   w₁ = x₁+y₁,   w₂ = x₂+y₂
        //   A  = x₁²+x₁,  B  = y₁²+y₁
        //
        //   x₃ = (d₁(x₁+x₂) + d₂·w₁·w₂ + A·(x₂(y₁+y₂+1)+y₁y₂))
        //         / (d₁ + A·w₂)
        //   y₃ = (d₁(y₁+y₂) + d₂·w₁·w₂ + B·(y₂(x₁+x₂+1)+x₁x₂))
        //         / (d₁ + B·w₂)
        // -----------------------------------------------------------------------

        fn add_binary(&self, other: &Self, curve: &EdwardsCurve<F>) -> Self {
            let w1 = &self.x + &self.y;
            let w2 = &other.x + &other.y;

            let x1_sq = <F as FieldOps>::square(&self.x);
            let a = &x1_sq + &self.x;

            let y1_sq = <F as FieldOps>::square(&self.y);
            let b = &y1_sq + &self.y;

            let w1w2 = &w1 * &w2;
            let d2_w1w2 = &curve.d2 * &w1w2;

            let x1x2 = &self.x * &other.x;
            let y1y2 = &self.y * &other.y;

            let x1_plus_x2 = &self.x + &other.x;
            let y1_plus_y2 = &self.y + &other.y;
            let one = F::one();

            let ysum_plus_one = &y1_plus_y2 + &one;
            let x_part = &other.x * &ysum_plus_one;
            let x_paren = &x_part + &y1y2;
            let ax = &a * &x_paren;
            let d1x = &curve.d1 * &x1_plus_x2;
            let x_num_tmp = &d1x + &d2_w1w2;
            let x_num = &x_num_tmp + &ax;
            let x_den = &curve.d1 + &(&a * &w2);

            let xsum_plus_one = &x1_plus_x2 + &one;
            let y_part = &other.y * &xsum_plus_one;
            let y_paren = &y_part + &x1x2;
            let by = &b * &y_paren;
            let d1y = &curve.d1 * &y1_plus_y2;
            let y_num_tmp = &d1y + &d2_w1w2;
            let y_num = &y_num_tmp + &by;
            let y_den = &curve.d1 + &(&b * &w2);

            let x_den_inv = <F as FieldOps>::invert(&x_den)
                .into_option()
                .expect("binary Edwards addition: x-denom must be invertible");
            let y_den_inv = <F as FieldOps>::invert(&y_den)
                .into_option()
                .expect("binary Edwards addition: y-denom must be invertible");

            let x3 = &x_num * &x_den_inv;
            let y3 = &y_num * &y_den_inv;

            Self::new(x3, y3)
        }

        // -----------------------------------------------------------------------
        // w-coordinate helpers for char-2 Montgomery ladder
        // -----------------------------------------------------------------------

        /// Differential addition on the `w`-line (`w = x + y`, char 2 only).
        ///
        /// Given `w₁ = w(Q−P)`, `w₂ = w(P)`, `w₃ = w(Q)`, compute `w₅ = w(P+Q)`.
        pub fn w_diff_add(w1: &F, w2: &F, w3: &F, curve: &EdwardsCurve<F>) -> F {
            let r = w2 * w3;
            let s = <F as FieldOps>::square(&r);

            let one = F::one();
            let tmp1 = &one + w2;
            let tmp2 = &tmp1 + w3;
            let rt = &r * &tmp2;
            let t = &rt + &s;

            let d1_inv = <F as FieldOps>::invert(&curve.d1)
                .into_option()
                .expect("d1 invertible");
            let d2_over_d1 = &curve.d2 * &d1_inv;
            let coeff = &d2_over_d1 + &one;

            let coeff_s = &coeff * &s;
            let den_tmp = &curve.d1 + &t;
            let den = &den_tmp + &coeff_s;
            let den_inv = <F as FieldOps>::invert(&den)
                .into_option()
                .expect("w-diff-add denominator invertible");

            let frac = &t * &den_inv;
            &frac + w1
        }

        /// `w`-coordinate doubling (`w = x + y`, char 2 only).
        ///
        /// Given `w₂ = w(P)`, compute `w₄ = w(2P)`.
        pub fn w_double(w2: &F, curve: &EdwardsCurve<F>) -> F {
            let a = <F as FieldOps>::square(w2);
            let j = <F as FieldOps>::square(&a);

            let d1_inv = <F as FieldOps>::invert(&curve.d1)
                .into_option()
                .expect("d1 invertible");
            let d2_over_d1 = &curve.d2 * &d1_inv;

            let num = &a + &j;
            let coeff_j = &d2_over_d1 * &j;
            let den_tmp = &curve.d1 + &a;
            let den = &den_tmp + &coeff_j;
            let den_inv = <F as FieldOps>::invert(&den)
                .into_option()
                .expect("w-double denominator invertible");

            &num * &den_inv
        }
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
