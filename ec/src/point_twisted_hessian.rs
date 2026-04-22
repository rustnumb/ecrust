//! Projective points on a twisted Hessian curve.
//!
//! We represent points on
//!
//! $$aX^3 + Y^3 + Z^3 = 3dXYZ$$
//!
//! by projective triples `(X:Y:Z)`.
//!
//! The formulas implemented here follow the twisted Hessian group law from
//! Decru--Kunzweiler (2026), §2.2 / Remark 2.11.

use core::fmt;

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::curve_twisted_hessian::TwistedHessianCurve;
use crate::point_ops::{PointAdd, PointOps};
use fp::field_ops::{FieldOps, FieldRandom};
use fp::{ref_field_impl, ref_field_trait_impl};

/// A projective point `(X:Y:Z)` on a twisted Hessian curve.
#[derive(Debug, Clone, Copy)]
pub struct TwistedHessianPoint<F: FieldOps> {
    /// Projective `X` coordinate.
    pub x: F,
    /// Projective `Y` coordinate.
    pub y: F,
    /// Projective `Z` coordinate.
    pub z: F,
}

impl<F> fmt::Display for TwistedHessianPoint<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_identity() {
            if f.alternate() {
                write!(f, "TwistedHessianPoint {{ O = (0:-1:1) }}")
            } else {
                write!(f, "O")
            }
        } else if self.is_zero_projective() {
            if f.alternate() {
                write!(f, "TwistedHessianPoint {{ invalid = (0:0:0) }}")
            } else {
                write!(f, "(0:0:0)")
            }
        } else if let Some((x_aff, y_aff)) = self.to_affine() {
            if f.alternate() {
                write!(
                    f,
                    "TwistedHessianPoint {{ X:Y:Z = ({}:{}:{}), x = {}, y = {} }}",
                    self.x, self.y, self.z, x_aff, y_aff
                )
            } else {
                write!(f, "({}, {})", x_aff, y_aff)
            }
        } else if f.alternate() {
            write!(
                f,
                "TwistedHessianPoint {{ X:Y:Z = ({}:{}:{}) }}",
                self.x, self.y, self.z
            )
        } else {
            write!(f, "({}:{}:{})", self.x, self.y, self.z)
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps> PartialEq for TwistedHessianPoint<F> {
        fn eq(&self, other: &Self) -> bool {
            let self_zero = self.is_zero_projective();
            let other_zero = other.is_zero_projective();
            if self_zero || other_zero {
                return self_zero && other_zero;
            }

            &self.x * &other.y == &other.x * &self.y
                && &self.x * &other.z == &other.x * &self.z
                && &self.y * &other.z == &other.y * &self.z
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps> Eq for TwistedHessianPoint<F> {}
}

impl<F: FieldOps> TwistedHessianPoint<F> {
    /// Construct a projective twisted Hessian point without validation.
    pub fn new(x: F, y: F, z: F) -> Self {
        Self { x, y, z }
    }

    /// Construct the affine point `(x, y)` as `(x:y:1)`.
    pub fn from_affine(x: F, y: F) -> Self {
        Self { x, y, z: F::one() }
    }

    /// Return the neutral element `(0:-1:1)`.
    pub fn identity() -> Self {
        let one = F::one();
        let minus_one = <F as FieldOps>::negate(&one);
        Self {
            x: F::zero(),
            y: minus_one,
            z: F::one(),
        }
    }

    /// Return `true` when the point is the neutral element.
    pub fn is_identity(&self) -> bool {
        if self.is_zero_projective() {
            return false;
        }

        bool::from(self.x.is_zero()) && F::add(&self.y, &self.z) == F::zero()
    }

    /// Return `true` if this is the invalid projective triple `(0:0:0)`.
    pub fn is_zero_projective(&self) -> bool {
        bool::from(self.x.is_zero())
            && bool::from(self.y.is_zero())
            && bool::from(self.z.is_zero())
    }

    /// Return `true` if the point lies on the line at infinity.
    pub fn is_at_infinity(&self) -> bool {
        bool::from(self.z.is_zero()) && !self.is_zero_projective()
    }

    /// Convert a finite projective point to affine coordinates.
    pub fn to_affine(&self) -> Option<(F, F)> {
        self.z.invert().into_option().map(|zinv| {
            (
                <F as FieldOps>::mul(&self.x, &zinv),
                <F as FieldOps>::mul(&self.y, &zinv),
            )
        })
    }
}

impl<F> ConditionallySelectable for TwistedHessianPoint<F>
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

ref_field_trait_impl! {
    impl<F: FieldOps + Copy + ConstantTimeEq> ConstantTimeEq for TwistedHessianPoint<F> {
        fn ct_eq(&self, other: &Self) -> Choice {
            let self_zero = self.is_zero_projective();
            let other_zero = other.is_zero_projective();
            if self_zero || other_zero {
                return Choice::from((self_zero && other_zero) as u8);
            }

            (&self.x * &other.y).ct_eq(&(&other.x * &self.y))
                & (&self.x * &other.z).ct_eq(&(&other.x * &self.z))
                & (&self.y * &other.z).ct_eq(&(&other.y * &self.z))
        }

        fn ct_ne(&self, other: &Self) -> Choice {
            !self.ct_eq(other)
        }
    }
}

ref_field_impl! {
    impl<F: FieldOps + FieldRandom> TwistedHessianPoint<F> {
        /// Negation on a twisted Hessian curve:
        ///
        /// $$-(X:Y:Z) = (X:Z:Y).$$
        pub fn negate(&self, _curve: &TwistedHessianCurve<F>) -> Self {
            Self::new(self.x.clone(), self.z.clone(), self.y.clone())
        }

        /// First addition formula from Remark 2.11.
        fn add_formula_1(&self, other: &Self) -> Self {
            let x1_sq = <F as FieldOps>::square(&self.x);
            let y1_sq = <F as FieldOps>::square(&self.y);
            let z1_sq = <F as FieldOps>::square(&self.z);
            let x2_sq = <F as FieldOps>::square(&other.x);
            let y2_sq = <F as FieldOps>::square(&other.y);
            let z2_sq = <F as FieldOps>::square(&other.z);

            let x3 = &(&x1_sq * &(&other.y * &other.z)) - &(&x2_sq * &(&self.y * &self.z));
            let y3 = &(&z1_sq * &(&other.x * &other.y)) - &(&z2_sq * &(&self.x * &self.y));
            let z3 = &(&y1_sq * &(&other.x * &other.z)) - &(&y2_sq * &(&self.x * &self.z));

            Self::new(x3, y3, z3)
        }

        /// Second addition formula from Remark 2.11.
        fn add_formula_2(&self, other: &Self, curve: &TwistedHessianCurve<F>) -> Self {
            let x1_sq = <F as FieldOps>::square(&self.x);
            let y1_sq = <F as FieldOps>::square(&self.y);
            let z1_sq = <F as FieldOps>::square(&self.z);
            let x2_sq = <F as FieldOps>::square(&other.x);
            let y2_sq = <F as FieldOps>::square(&other.y);
            let z2_sq = <F as FieldOps>::square(&other.z);

            let x3 = &(&z2_sq * &(&self.x * &self.z)) - &(&y1_sq * &(&other.x * &other.y));
            let y3 = &(&y2_sq * &(&self.y * &self.z)) - &(&curve.a * &(&x1_sq * &(&other.x * &other.z)));
            let z3 = &(&curve.a * &(&x2_sq * &(&self.x * &self.y))) - &(&z1_sq * &(&other.y * &other.z));

            Self::new(x3, y3, z3)
        }

        /// Add two projective twisted Hessian points.
        ///
        /// The first branch is used generically; the second branch covers the
        /// complementary exceptional set.
        pub fn add(&self, other: &Self, curve: &TwistedHessianCurve<F>) -> Self {
            if self.is_identity() {
                return other.clone();
            }
            if other.is_identity() {
                return self.clone();
            }

            let r = self.add_formula_1(other);
            if !r.is_zero_projective() {
                return r;
            }

            let s = self.add_formula_2(other, curve);
            if !s.is_zero_projective() {
                return s;
            }

            if other == &self.negate(curve) {
                return Self::identity();
            }

            panic!("twisted Hessian addition failed for valid-looking inputs; both unified formula branches vanished");
        }

        /// Doubling via the unified addition formulas.
        pub fn double(&self, curve: &TwistedHessianCurve<F>) -> Self {
            self.add(self, curve)
        }

        /// Variable-time double-and-add scalar multiplication.
        pub fn scalar_mul(&self, k: &[u64], curve: &TwistedHessianCurve<F>) -> Self {
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
    impl<F: FieldOps> PointOps for TwistedHessianPoint<F> {
        type BaseField = F;
        type Curve = TwistedHessianCurve<F>;

        fn identity(_curve: &Self::Curve) -> Self {
            TwistedHessianPoint::<F>::identity()
        }

        fn is_identity(&self) -> bool {
            TwistedHessianPoint::<F>::is_identity(self)
        }

        fn negate(&self, curve: &Self::Curve) -> Self {
            TwistedHessianPoint::<F>::negate(self, curve)
        }

        fn scalar_mul(&self, k: &[u64], curve: &Self::Curve) -> Self {
            TwistedHessianPoint::<F>::scalar_mul(self, k, curve)
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps> PointAdd for TwistedHessianPoint<F> {
        fn add(&self, other: &Self, curve: &Self::Curve) -> Self {
            TwistedHessianPoint::<F>::add(self, other, curve)
        }
    }
}
