//! Elliptic curve definition in Edwards form.
//!
//! # Odd characteristic
//!
//! ```text
//! x² + y² = 1 + d·x²·y²
//! ```
//!
//! Parameters: `d1 = 0` (unused), `d2 = d`.
//! Identity: `(0, 1)`.  Complete when `d` is not a square.
//!
//! # Characteristic 2
//!
//! ```text
//! d₁(x + y) + d₂(x² + y²) = xy + xy(x + y) + x²y²
//! ```
//!
//! Parameters: `d1`, `d2` with `d₁ ≠ 0`, `d₂ ≠ d₁² + d₁`.
//! Identity: `(0, 0)`.  Complete when `Tr(d₂) = 1`.
//!
//! Reference (char 2): Bernstein–Lange–Rezaeian Farashahi,
//!   "Binary Edwards Curves", 2008.
//! Reference (odd char): <https://hyperelliptic.org/EFD/g1p/auto-edwards.html>

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::field_ops::FieldOps;
use fp::{ref_field_impl, ref_field_trait_impl};

use crate::curve_ops::Curve;
use crate::point_edwards::EdwardsPoint;

/// An Edwards curve over a field `F`, covering both odd and even characteristic.
///
/// In odd characteristic only `d2` is used (the parameter `d`).
/// In characteristic 2 both `d1` and `d2` are used.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct EdwardsCurve<F: FieldOps> {
    pub d1: F,
    pub d2: F,
}

ref_field_impl! {
    impl<F> EdwardsCurve<F> {
        /// Construct an odd-characteristic Edwards curve `x² + y² = 1 + d·x²·y²`.
        ///
        /// Stores `d` as `d2`; `d1` is set to zero (unused).
        pub fn new(d: F) -> Self {
            assert!(F::characteristic()[0] != 2, "use new_binary() for char 2");
            assert!(d != F::zero(), "d must be nonzero");
            assert!(d != F::one(), "d must not be 1");
            Self { d1: F::zero(), d2: d }
        }

        /// Construct a binary Edwards curve
        /// `d₁(x+y) + d₂(x²+y²) = xy + xy(x+y) + x²y²`.
        pub fn new_binary(d1: F, d2: F) -> Self {
            assert_eq!(F::characteristic()[0], 2, "use new() for odd char");
            assert!(d1 != F::zero(), "d1 must be nonzero");

            let d1_sq = <F as FieldOps>::square(&d1);
            let forbidden = &d1_sq + &d1;
            assert!(d2 != forbidden, "d2 must differ from d1² + d1");

            Self { d1, d2 }
        }

        /// Convenience accessor: the parameter `d` in odd characteristic.
        pub fn d(&self) -> F {
            self.d2.clone()
        }

        /// Check whether the affine point `(x, y)` lies on the curve.
        pub fn contains(&self, x: &F, y: &F) -> bool {
            if F::characteristic()[0] != 2 {
                // x² + y² == 1 + d·x²·y²
                let x2 = <F as FieldOps>::square(x);
                let y2 = <F as FieldOps>::square(y);

                let lhs = &x2 + &y2;

                let x2y2 = &x2 * &y2;
                let dx2y2 = &self.d2 * &x2y2;
                let one = F::one();
                let rhs = &one + &dx2y2;

                lhs == rhs
            } else {
                // d₁(x+y) + d₂(x²+y²) == xy + xy(x+y) + x²y²
                let x2 = <F as FieldOps>::square(x);
                let y2 = <F as FieldOps>::square(y);
                let xy = x * y;
                let xpy = x + y;

                let lhs1 = &self.d1 * &xpy;
                let x2py2 = &x2 + &y2;
                let lhs2 = &self.d2 * &x2py2;
                let lhs = &lhs1 + &lhs2;

                let xy_xpy = &xy * &xpy;
                let x2y2 = &x2 * &y2;
                let rhs_tmp = &xy + &xy_xpy;
                let rhs = &rhs_tmp + &x2y2;

                lhs == rhs
            }
        }
    }
}

ref_field_trait_impl! {
    impl<F> Curve for EdwardsCurve<F> {
        type BaseField = F;
        type Point = EdwardsPoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            self.contains(&point.x, &point.y)
        }

        fn random_point(&self) -> Self::Point {
            todo!()
        }

        fn j_invariant(&self) -> F {
            if F::characteristic()[0] != 2 {
                // j = 16(1 + 14d + d²)³ / (d(1 − d)⁴)
                let d_sq = <F as FieldOps>::square(&self.d2);

                let one = F::one();
                let two = <F as FieldOps>::double(&one);
                let four = <F as FieldOps>::double(&two);
                let eight = <F as FieldOps>::double(&four);
                let fourteen_tmp = &eight + &four;
                let fourteen = &fourteen_tmp + &two;
                let sixteen = <F as FieldOps>::double(&eight);

                let fourteen_d = &fourteen * &self.d2;
                let inner_tmp = &one + &fourteen_d;
                let inner = &inner_tmp + &d_sq;
                let inner_sq = <F as FieldOps>::square(&inner);
                let inner_cubed = &inner * &inner_sq;
                let numer = &sixteen * &inner_cubed;

                let one_minus_d = &one - &self.d2;
                let omd2 = <F as FieldOps>::square(&one_minus_d);
                let omd4 = <F as FieldOps>::square(&omd2);
                let denom = &self.d2 * &omd4;

                let denom_inv = <F as FieldOps>::invert(&denom)
                    .into_option()
                    .expect("d(1-d)^4 must be invertible");

                &numer * &denom_inv
            } else {
                // j = 1 / (d₁⁴ (d₁⁴ + d₁² + d₂²))
                let d1_sq = <F as FieldOps>::square(&self.d1);
                let d1_4 = <F as FieldOps>::square(&d1_sq);
                let d2_sq = <F as FieldOps>::square(&self.d2);

                let inner_tmp = &d1_4 + &d1_sq;
                let inner = &inner_tmp + &d2_sq;
                let denom = &d1_4 * &inner;

                <F as FieldOps>::invert(&denom)
                    .into_option()
                    .expect("j-invariant denominator must be invertible")
            }
        }

        fn a_invariants(&self) -> Vec<Self::BaseField> {
            if F::characteristic()[0] != 2 {
                vec![self.d2.clone()]
            } else {
                vec![self.d1.clone(), self.d2.clone()]
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for EdwardsCurve<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            d1: F::conditional_select(&a.d1, &b.d1, choice),
            d2: F::conditional_select(&a.d2, &b.d2, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.d1, &other.d1, choice);
        F::conditional_assign(&mut self.d2, &other.d2, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.d1, &mut b.d1, choice);
        F::conditional_swap(&mut a.d2, &mut b.d2, choice);
    }
}

impl<F> ConstantTimeEq for EdwardsCurve<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.d1.ct_eq(&other.d1) & self.d2.ct_eq(&other.d2)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}