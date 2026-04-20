//! Point representation and group law for Legendre curves.
//!
//! We work with the Legendre family
//!
//! $$
//! E_\lambda : y^2 = x(x-1)(x-\lambda)
//! $$
//!
//! over fields of odd characteristic.
//!
//! # Representation
//!
//! Points are stored in affine coordinates `(x, y)` together with a distinguished
//! point at infinity `O`, represented by `infinity = true`.
//!
//! # Group law
//!
//! Writing the curve as
//!
//! $$
//! y^2 = x^3 - (1+\lambda)x^2 + \lambda x,
//! $$
//!
//! the affine group law is the usual chord-and-tangent law for a Weierstrass
//! model with an `x²` term.
//!
//! - Identity: `O`
//! - Negation: `-(x,y) = (x,-y)`
//! - Visible 2-torsion:
//!
//! $$
//! (0,0),\quad (1,0),\quad (\lambda,0).
//! $$
//!
//! For distinct finite points `P=(x₁,y₁)` and `Q=(x₂,y₂)` with `x₁ != x₂`:
//!
//! $$
//! m = \frac{y_2-y_1}{x_2-x_1},
//! $$
//!
//! $$
//! x_3 = m^2 + (1+\lambda) - x_1 - x_2,
//! \qquad
//! y_3 = m(x_1-x_3)-y_1.
//! $$
//!
//! For doubling `P=(x₁,y₁)` with `y₁ != 0`:
//!
//! $$
//! m = \frac{3x_1^2 - 2(1+\lambda)x_1 + \lambda}{2y_1},
//! $$
//!
//! $$
//! x_3 = m^2 + (1+\lambda) - 2x_1,
//! \qquad
//! y_3 = m(x_1-x_3)-y_1.
//! $$
//!
//! # References
//!
//! - Roland Auer and Jaap Top, *Legendre Elliptic Curves over Finite Fields*,
//!   Journal of Number Theory 95 (2002), 303–312.
//! - HongFeng Wu and RongQuan Feng,
//!   *On the isomorphism classes of Legendre elliptic curves over finite fields*,
//!   Sci. China Math. 54(9) (2011), 1885–1890.

use core::fmt;
//use std::os::unix::raw::ino_t;
//use crypto_bigint::modular::ConstMontyForm;
use crate::curve_legendre::LegendreCurve;
use crate::point_ops::PointOps;
use fp::field_ops::FieldOps;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::ref_field_impl;

/// A point on a Legendre elliptic curve over `F`.
///
/// The identity element is the point at infinity. It is represented by
/// `infinity = true`; in that case the `x` and `y` fields are dummy values.
#[derive(Debug, Clone, Copy)]
pub struct LegendrePoint<F: FieldOps> {
    /// Affine x-coordinate.
    pub x: F,
    /// Affine y-coordinate.
    pub y: F,
    /// Whether this point is the point at infinity.
    pub infinity: bool,
}

impl<F: FieldOps> PartialEq for LegendrePoint<F>
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

impl<F> fmt::Display for LegendrePoint<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.infinity {
            if f.alternate() {
                write!(f, "LegendrePoint {{ O }}")
            } else {
                write!(f, "O")
            }
        } else if f.alternate() {
            write!(f, "LegendrePoint {{ x = {}, y = {} }}", self.x, self.y)
        } else {
            write!(f, "({}, {})", self.x, self.y)
        }
    }
}

impl<F: FieldOps> LegendrePoint<F> {
    /// Constructs a finite affine point `(x, y)`.
    ///
    /// No on-curve check is performed. Use
    /// [`LegendreCurve::contains`] or [`crate::curve_ops::Curve::is_on_curve`]
    /// if validation is needed.
    pub fn new(x: F, y: F) -> Self {
        Self { x, y, infinity: false }
    }

    /// Returns the identity element `O`.
    ///
    /// Internally this is represented by `infinity = true`; the stored
    /// coordinates are dummy values.
    pub fn identity() -> Self {
        Self {
            x: F::zero(),
            y: F::one(),
            infinity: true,
        }
    }

    /// Returns `true` if this point is the identity.
    pub fn is_identity(&self) -> bool {
        self.infinity
    }
}

impl<F> ConditionallySelectable for LegendrePoint<F>
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

impl<F> ConstantTimeEq for LegendrePoint<F>
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

ref_field_impl!{
    impl<F: FieldOps> LegendrePoint<F> {
        ///
        /// For a finite affine point `(x, y)` on a Legendre curve, the inverse is
        /// `(x, -y)`. The point at infinity is its own inverse.
        pub fn negate(&self, _curve: &LegendreCurve<F>) -> Self {
            if self.infinity {
                return Self::identity();
            }
            Self::new(self.x, -self.y)
        }

        /// Doubles the point: returns `[2]P`.
        ///
        /// If `P = O`, returns `O`.
        /// If `y = 0`, then `P` is a nontrivial 2-torsion point and `[2]P = O`.
        pub fn double(&self, curve: &LegendreCurve<F>) -> Self {
            if self.infinity {
                return Self::identity();
            }

            let denom = <F as FieldOps>::double(&self.y);
            let denom_inv = match denom.invert().into_option() {
                Some(inv) => inv,
                None => return Self::identity(),
            };

            let x1_sq = <F as FieldOps>::square(&self.x);
            let three_x1_sq = &F::from_u64(3) * &x1_sq;
            let one_plus_lambda = &F::one() + &curve.lambda;
            let twice_one_plus_lambda = &F::from_u64(2) * &one_plus_lambda;
            let num = &(&three_x1_sq - &(&twice_one_plus_lambda * &self.x)) + &curve.lambda;

            let m = &num * &denom_inv;

            let x3 = &<F as FieldOps>::square(&m) + &(&one_plus_lambda - &(&F::from_u64(2) * &self.x));
            let y3 = &(&m * &(&self.x - &x3)) - &self.y;

            Self::new(x3, y3)
        }

        /// Adds two points: returns `P + Q`.
        ///
        /// Handles the usual special cases:
        /// - `O + Q = Q`
        /// - `P + O = P`
        /// - `P = Q` uses doubling
        /// - `P = -Q` returns `O`
        ///
        /// For distinct finite points `P = (x1, y1)` and `Q = (x2, y2)` with
        /// `x1 != x2`, the slope is
        ///
        /// `m = (y2 - y1) / (x2 - x1)`
        ///
        /// and
        ///
        /// `x3 = m^2 + (1+λ) - x1 - x2`
        ///
        /// `y3 = m(x1 - x3) - y1`.
        pub fn add(&self, other: &Self, curve: &LegendreCurve<F>) -> Self {
            if self.infinity {
                return *other;
            }
            if other.infinity {
                return *self;
            }

            if self.x == other.x {
                if self.y == other.y {
                    return self.double(curve);
                }
                return Self::identity();
            }

            let dx = &other.x - &self.x;
            let dy = &other.y - &self.y;

            let dx_inv = dx
                .invert()
                .into_option()
                .expect("x1 != x2, so dx must be invertible");

            let m = &dy * &dx_inv;
            let one_plus_lambda = &F::one() + &curve.lambda;

            let x3 = &<F as FieldOps>::square(&m) + &(&one_plus_lambda - &(&self.x - &other.x));
            let y3 = &(&m * &(&self.x - &x3)) - &self.y;

            Self::new(x3, y3)
        }

        /// Multiply `self` by `k`
        ///
        /// # Arguments
        ///
        /// * `&self` - Point on curve (type: `Self`)
        /// * `k` - Integer (type: `&[u64]`)
        /// * `curve` - The curve we're on (type: `&<LegendrePoint<F> as PointOps>::Curve`)
        ///
        /// # Returns
        ///
        /// The point `k * self` (type: `Self`)
        pub fn scalar_mul(&self, k: &[u64], curve: &<LegendrePoint<F> as PointOps>::Curve) -> Self {
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

impl<F> PointOps for LegendrePoint<F>
where
    F: FieldOps,
{
    type BaseField = F;
    type Curve = LegendreCurve<F>;

    fn identity(_curve: &Self::Curve) -> Self {
        LegendrePoint::<F>::identity()
    }

    fn is_identity(&self) -> bool {
        self.infinity
    }

    fn negate(&self, curve: &Self::Curve) -> Self {
        LegendrePoint::<F>::negate(self, curve)
    }

    fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
        LegendrePoint::<F>::scalar_mul(self, k, curve)
    }
}

impl<F> crate::point_ops::PointAdd for LegendrePoint<F>
where
    F: FieldOps,
{
    fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
        LegendrePoint::<F>::add(self, other, curve)
    }
}
