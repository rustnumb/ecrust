//! Affine points on a Jacobi quartic curve.
//!
//! We use the affine model
//!
//! $$
//! y^2 = d x^4 + 2 a x^2 + 1
//! $$
//!
//! with neutral element $(0, 1)$ and negation $-(x, y) = (-x, y)$.
//! The doubling formulas are taken from equations (9) and (10) in
//! *Jacobi Quartic Curves Revisited*; the general addition uses affine
//! formulas (1) and (2).
//!
//! Important: these are affine formulas. Like the existing Edwards code in this
//! crate, they assume denominators are invertible. For exceptional inputs
//! where the result leaves this affine chart, the code panics instead of trying
//! to model the desingularized points at infinity.
use core::fmt;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_jacobi_quartic::JacobiQuarticCurve;
use crate::point_ops::{PointAdd, PointOps};
use fp::field_ops::FieldOps;
use fp::{ref_field_impl, ref_field_trait_impl};

/// An affine point `(x, y)` on a Jacobi quartic curve.
#[derive(Debug, Clone, Copy)]
pub struct JacobiQuarticPoint<F: FieldOps> {
    /// The x coordinate of a point
    pub x: F,
    /// The y coordinate of a point
    pub y: F,
}

impl<F> fmt::Display for JacobiQuarticPoint<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_identity() {
            if f.alternate() {
                write!(f, "JacobiQuarticPoint {{ O = (0,1) }}")
            } else {
                write!(f, "O")
            }
        } else if f.alternate() {
            write!(f, "JacobiQuarticPoint {{ x = {}, y = {} }}", self.x, self.y)
        } else {
            write!(f, "({}, {})", self.x, self.y)
        }
    }
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
        Self {
            x: F::zero(),
            y: F::one(),
        }
    }

    pub fn is_identity(&self) -> bool {
        self.x == F::zero() && self.y == F::one()
    }

    /// The affine order-2 point `(0, -1)`.
    pub fn order_two_point() -> Self {
        let one = F::one();
        let minus_one = <F as FieldOps>::negate(&one);
        Self {
            x: F::zero(),
            y: minus_one,
        }
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

ref_field_impl! {
    impl<F> JacobiQuarticPoint<F> {
        /// Negation on a Jacobi quartic: `-(x, y) = (-x, y)`.
        pub fn negate(&self, _curve: &JacobiQuarticCurve<F>) -> Self {
            Self::new(-&self.x, self.y.clone())
        }

        /// Dedicated affine doubling from equations (9) and (10):
        ///
        /// ```text
        /// μ  = 2y / (2 + 2ax² − y²)
        /// x₃ = μx
        /// y₃ = μ(μ − y) − 1.
        /// ```
        pub fn double(&self, curve: &JacobiQuarticCurve<F>) -> Self {
            if self.is_identity() {
                return *self;
            }

            let x2 = <F as FieldOps>::square(&self.x);
            let y2 = <F as FieldOps>::square(&self.y);

            let one = F::one();
            let two = <F as FieldOps>::double(&one);

            let ax2 = &curve.a * &x2;
            let two_ax2 = &two * &ax2;
            let denom_tmp = &two + &two_ax2;
            let denom = &denom_tmp - &y2;
            let denom_inv = <F as FieldOps>::invert(&denom)
                .into_option()
                .expect("Jacobi quartic doubling denominator vanished; result leaves affine chart or input is exceptional");

            let two_y = &two * &self.y;
            let mu = &two_y * &denom_inv;
            let x3 = &mu * &self.x;

            let mu_minus_y = &mu - &self.y;
            let mu_times = &mu * &mu_minus_y;
            let y3 = &mu_times - &one;

            Self::new(x3, y3)
        }

        /// Affine addition using equations (1) and (2):
        ///
        /// ```text
        /// x₃ = (x₁y₂ + y₁x₂)/(1 − d x₁²x₂²)
        /// y₃ = ((y₁y₂ + 2ax₁x₂)(1 + d x₁²x₂²) + 2d x₁x₂(x₁² + x₂²))
        ///      /(1 − d x₁²x₂²)².
        /// ```
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

            let x1_sq = <F as FieldOps>::square(&self.x);
            let x2_sq = <F as FieldOps>::square(&other.x);
            let x1x2 = &self.x * &other.x;
            let y1y2 = &self.y * &other.y;
            let x1_sq_x2_sq = &x1_sq * &x2_sq;
            let dx4 = &curve.d * &x1_sq_x2_sq;

            let one = F::one();
            let two = <F as FieldOps>::double(&one);

            let denom = &one - &dx4;
            let denom_inv = <F as FieldOps>::invert(&denom)
                .into_option()
                .expect("Jacobi quartic addition denominator vanished; choose d nonsquare / odd-order subgroup or use a projective model");

            let x1y2 = &self.x * &other.y;
            let y1x2 = &self.y * &other.x;
            let x_num = &x1y2 + &y1x2;
            let x3 = &x_num * &denom_inv;

            let two_a = &two * &curve.a;
            let two_a_x1x2 = &two_a * &x1x2;
            let y_part = &y1y2 + &two_a_x1x2;
            let one_plus_dx4 = &one + &dx4;
            let first_term = &y_part * &one_plus_dx4;

            let x_sq_sum = &x1_sq + &x2_sq;
            let two_d = &two * &curve.d;
            let two_d_x1x2 = &two_d * &x1x2;
            let second_term = &two_d_x1x2 * &x_sq_sum;

            let numer_y = &first_term + &second_term;
            let denom_inv_sq = <F as FieldOps>::square(&denom_inv);
            let y3 = &numer_y * &denom_inv_sq;

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
}

ref_field_trait_impl! {
    impl<F> PointOps for JacobiQuarticPoint<F> {
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
}

ref_field_trait_impl! {
    impl<F> PointAdd for JacobiQuarticPoint<F> {
        fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
            JacobiQuarticPoint::<F>::add(self, other, curve)
        }
    }
}