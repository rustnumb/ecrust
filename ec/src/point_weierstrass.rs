//! Elliptic curve point representation and group law.
//!
//! # Representation
//!
//! Points are stored in **affine** coordinates $(x, y)$ or as the
//! distinguished point at infinity $O$ (the group identity).
//!
//! # Group law
//!
//! The addition formulas implement the general Weierstrass group law:
//!
//! | Operation         | Algorithm                                | Cost           |
//! |-------------------|------------------------------------------|----------------|
//! | Negate            | `-(x,y) = (x, -y - a_1x - a_3)`          | 3 mul + 2 add  |
//! | Add  (P ≠ ±Q)     | Chord-and-tangent (general Weierstrass)  | 1 inv + 6 mul  |
//! | Double            | Tangent (general Weierstrass)            | 1 inv + 7 mul  |
//! | Scalar multiply   | Montgomery ladder (scalar-constant-time) | O(n) doubles   |
//! | Operation        | Algorithm                                | Cost              |
//! |------------------|------------------------------------------|-------------------|
//! | Negate           | $-(x,y) = (x, -y - a_1 x - a_3)$         | $3$ mul + $2$ add |
//! | Add $(P \ne \pm Q)$ | Chord-and-tangent (general Weierstrass) | $1$ inv + $6$ mul |
//! | Double           | Tangent (general Weierstrass)            | $1$ inv + $7$ mul |
//! | Scalar multiply  | Montgomery ladder (constant-time scalar) | $O(n)$ doublings  |

// WARNING: SOME OF THE FUNCTIONS BELOW USE BRANCHES DEPENDING
// ON WHETHER A POINT IS AT INFINITY OR NOT!!!!

use core::fmt;
//use std::os::unix::raw::ino_t;
//use crypto_bigint::modular::ConstMontyForm;
use crate::curve_weierstrass::WeierstrassCurve;
use crate::point_ops::PointOps;
use fp::field_ops::FieldOps;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::{ref_field_impl, ref_field_trait_impl, ref_field_trait_impl_path};

/// An affine point on a Weierstrass elliptic curve over `F`.
///
/// The point at infinity is represented by `infinity = true`; in that case
/// the `x` and `y` fields are meaningless (set to zero by convention).
#[derive(Debug, Clone, Copy)]
pub struct AffinePoint<F: FieldOps> {
    /// x coordinate
    pub x: F,
    /// y coordinate
    pub y: F,
    /// true if and only if the ponit at infinity
    pub infinity: bool,
}

// ---------------------------------------------------------------------------
// Manual trait impls  (avoid over-constraining with #[derive])
// ---------------------------------------------------------------------------

impl<F: FieldOps> PartialEq for AffinePoint<F>
where
    F: FieldOps + ConstantTimeEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self.infinity, other.infinity) {
            (true, true) => true,
            (true, false) => false,
            (false, true) => false,
            (false, false) => self.x == other.x && self.y == other.y,
        }
    }
}

impl<F> fmt::Display for AffinePoint<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.infinity {
            if f.alternate() {
                write!(f, "AffinePoint {{ O }}")
            } else {
                write!(f, "O")
            }
        } else if f.alternate() {
            write!(f, "AffinePoint {{ x = {}, y = {} }}", self.x, self.y)
        } else {
            write!(f, "({}, {})", self.x, self.y)
        }
    }
}



impl<F: FieldOps> Eq for AffinePoint<F>
where
    F: FieldOps + ConstantTimeEq,
{
}

impl<F> core::hash::Hash for AffinePoint<F>
where
    F: FieldOps + ConstantTimeEq + core::hash::Hash,
{
    //  Must agree with `PartialEq`: two identities compare equal regardless
    //  of their (unused) `x`, `y` fields, so they must hash the same.
    //  We fold `infinity` in first and only hash coordinates for finite
    //  points.
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.infinity.hash(state);
        if !self.infinity {
            self.x.hash(state);
            self.y.hash(state);
        }
    }
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<F: FieldOps> AffinePoint<F> {
    /// Construct a finite affine point `(x, y)`.
    ///
    /// **No** on-curve check is performed; use
    /// [`WeierstrassCurve::contains`] if you need validation.
    pub fn new(x: F, y: F) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }

    /// The point at infinity `O` (the group identity).
    pub fn identity() -> Self {
        Self {
            x: F::zero(),
            y: F::zero(),
            infinity: true,
        }
    }

    /// Returns `true` if this is the point at infinity.
    pub fn is_identity(&self) -> bool {
        self.infinity
    }
}

// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for AffinePoint<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let ai = a.infinity as u8;
        let bi = b.infinity as u8;
        let infinity = u8::conditional_select(&ai, &bi, choice) != 0;

        Self {
            x: F::conditional_select(&a.x, &b.x, choice),
            y: F::conditional_select(&a.y, &b.y, choice),
            infinity,
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.x, &other.x, choice);
        F::conditional_assign(&mut self.y, &other.y, choice);

        let mut inf = self.infinity as u8;
        let other_inf = other.infinity as u8;
        inf.conditional_assign(&other_inf, choice);
        self.infinity = inf != 0;
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.x, &mut b.x, choice);
        F::conditional_swap(&mut a.y, &mut b.y, choice);

        let mut ai = a.infinity as u8;
        let mut bi = b.infinity as u8;
        u8::conditional_swap(&mut ai, &mut bi, choice);
        a.infinity = ai != 0;
        b.infinity = bi != 0;
    }
}

impl<F> ConstantTimeEq for AffinePoint<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        let self_inf = Choice::from(self.infinity as u8);
        let other_inf = Choice::from(other.infinity as u8);

        let both_inf = self_inf & other_inf;
        let both_finite = !self_inf & !other_inf;
        let coords_eq = self.x.ct_eq(&other.x) & self.y.ct_eq(&other.y);

        both_inf | (both_finite & coords_eq)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}

// ---------------------------------------------------------------------------
// Group operations
// ---------------------------------------------------------------------------

ref_field_impl! {
    impl<F> AffinePoint<F> {
        /// Negate a point:  `-(x, y) = (x, −y − a₁x − a₃)`.
        pub fn negate(&self, curve: &WeierstrassCurve<F>) -> Self {
            if self.infinity {
                return Self::identity();
            }

            let neg_y = -&(&self.y + &(&(&curve.a1 * &self.x) - &curve.a3));

            Self::new(self.x.clone(), neg_y)
        }

        /// Double a point:  `[2]P`.
        ///
        /// Uses the tangent-line formula for the general Weierstrass model:
        ///
        /// ```text
        /// λ = (3x₁² + 2a₂x₁ + a₄ − a₁y₁) / (2y₁ + a₁x₁ + a₃)
        /// x₃ = λ² + a₁λ − a₂ − 2x₁
        /// y₃ = λ(x₁ − x₃) − y₁ − a₁x₃ − a₃
        /// ```
        ///
        /// Returns `O` when the tangent is vertical (i.e. `2y + a₁x + a₃ = 0`).
        pub fn double(&self, curve: &WeierstrassCurve<F>) -> Self {
        if self.infinity {
            return Self::identity();
        }

        // denominator = 2y₁ + a₁x₁ + a₃
        let two_y = <F as FieldOps>::double(&self.y);
        let a1x = &curve.a1 * &self.x;
        let tmp = &two_y + &a1x;
        let denom = &tmp + &curve.a3;

        // If denominator is zero the tangent is vertical → result is O.
        let denom_inv = match denom.invert().into_option() {
            Some(inv) => inv,
            None => return Self::identity(),
        };

        // numerator = 3x₁² + 2a₂x₁ + a₄ − a₁y₁
        let numer = {
                let x1_sq = <F as FieldOps>::square(&self.x);
                let two_x1_sq = <F as FieldOps>::double(&x1_sq);
                let three_x1_sq = &x1_sq + &two_x1_sq;

                let a2x1 = &curve.a2 * &self.x;
                let two_a2_x1 = <F as FieldOps>::double(&a2x1);

                let a1y1 = &curve.a1 * &self.y;

                let tmp1 = &three_x1_sq + &two_a2_x1;
                let tmp2 = &tmp1 + &curve.a4;
                &tmp2 - &a1y1
            };

            let lambda = &numer * &denom_inv;

            // x₃ = λ² + a₁λ − a₂ − 2x₁
            let x3 = {
                let lam_sq = <F as FieldOps>::square(&lambda);
                let a1_lam = &curve.a1 * &lambda;
                let two_x1 = <F as FieldOps>::double(&self.x);

                let tmp1 = &lam_sq + &a1_lam;
                let tmp2 = &tmp1 - &curve.a2;
                &tmp2 - &two_x1
            };

            // y₃ = λ(x₁ − x₃) − y₁ − a₁x₃ − a₃
            let y3 = {
                let dx = &self.x - &x3;
                let lam_dx = &lambda * &dx;
                let a1x3 = &curve.a1 * &x3;

                let tmp1 = &lam_dx - &self.y;
                let tmp2 = &tmp1 - &a1x3;
                &tmp2 - &curve.a3
            };

            Self::new(x3, y3)
        }

        /// Add two points:  `P + Q`.
        ///
        /// Handles all cases:
        ///   - Either operand is `O` → return the other.
        ///   - `P = Q` → delegate to [`double`](Self::double).
        ///   - `P = −Q` (same x, opposite y) → return `O`.
        ///   - General chord:
        ///
        /// ```text
        /// λ  = (y₂ − y₁) / (x₂ − x₁)
        /// x₃ = λ² + a₁λ − a₂ − x₁ − x₂
        /// y₃ = λ(x₁ − x₃) − y₁ − a₁x₃ − a₃
        /// ```
        pub fn add(&self, other: &Self, curve: &WeierstrassCurve<F>) -> Self {
            // O + Q = Q
            if self.infinity {
                return other.clone();
            }
            // P + O = P
            if other.infinity {
                return self.clone();
            }

            // Same x-coordinate?
            if self.x == other.x {
                if self.y == other.y {
                    // P = Q  → doubling
                    return self.double(curve);
                }
                // P = −Q  → identity
                // (This also covers the char-2 case where y₁ + y₂ + a₁x + a₃ = 0.)
                return Self::identity();
            }

            // General chord
            let dx = &other.x - &self.x;
            let dy = &other.y - &self.y;

            // dx ≠ 0 guaranteed by the x₁ ≠ x₂ check above
            let dx_inv = <F as FieldOps>::invert(&dx)
                .into_option()
                .expect("dx must be invertible (x₁ ≠ x₂)");
            let lambda = &dy * &dx_inv;

            // x₃ = λ² + a₁λ − a₂ − x₁ − x₂
            let x3 = {
                let lam_sq = <F as FieldOps>::square(&lambda);
                let a1_lam = &curve.a1 * &lambda;

                let tmp1 = &lam_sq + &a1_lam;
                let tmp2 = &tmp1 - &curve.a2;
                let tmp3 = &tmp2 - &self.x;
                &tmp3 - &other.x
            };

            // y₃ = λ(x₁ − x₃) − y₁ − a₁x₃ − a₃
            let y3 = {
                let dx = &self.x - &x3;
                let lam_dx = &lambda * &dx;
                let a1x3 = &curve.a1 * &x3;

                let tmp1 = &lam_dx - &self.y;
                let tmp2 = &tmp1 - &a1x3;
                &tmp2 - &curve.a3
            };

            Self::new(x3, y3)

        }

        /// Multiply `self` by `k`
        ///
        /// # Arguments
        ///
        /// * `&self` - Point on curve (type: `Self`)
        /// * `k` - Integer (type: `&[u64]`)
        /// * `curve` - The curve we're on (type: `&<AffinePoint<F> as PointOps>::Curve`)
        ///
        /// # Returns
        ///
        /// The point `k * self` (type: `Self`)
        pub fn scalar_mul(&self, k: &[u64], curve: &<AffinePoint<F> as PointOps>::Curve) -> Self {
            let mut r0 = Self::identity();
            let mut r1 = self.clone();

            for &limb in k.iter().rev() {
                for bit in (0..64).rev() {
                    let choice = Choice::from(((limb >> bit) & 1) as u8);

                    Self::conditional_swap(&mut r0, &mut r1, choice);

                    let sum = r0.add(&r1, curve);
                    let dbl = r0.double(curve);
                    r1 = sum;
                    r0 = dbl;

                    Self::conditional_swap(&mut r0, &mut r1, choice);
                }
            }

            r0
        }
    }
}

ref_field_trait_impl! {
    impl<F> PointOps for AffinePoint<F> {
        type BaseField = F;
        type Curve = WeierstrassCurve<F>;

        fn identity(_curve: &Self::Curve) -> Self {
            AffinePoint::<F>::identity()
        }

        fn is_identity(&self) -> bool {
            self.infinity
        }

        fn negate(&self, curve: &Self::Curve) -> Self {
            AffinePoint::<F>::negate(self, curve)
        }

        fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
            AffinePoint::<F>::scalar_mul(self, k, curve)
        }
    }
}

ref_field_trait_impl_path! {
    impl<F> (crate::point_ops::PointAdd) for AffinePoint<F> {
        fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
            AffinePoint::<F>::add(self, other, curve)
        }
    }
}