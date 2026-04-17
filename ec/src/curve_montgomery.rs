//! Elliptic curve definition (Montgomery form).
//!
//! # Equation
//!
//! Given a field $F$ of odd characteristic, a Montgomery curve over $F$ is given by
//!
//! $$
//! B y^2 = x(x^2 + A x + 1)
//! $$
//!
//! where $B \ne 0$.
//!
//! Given a binary field $F$, a Montgomery curve over $F$ is given by
//!
//! $$
//! y^2 + xy = x(x^2 + A x + B^2)
//! $$
//!
//! where $B \ne 0$.
//!
//! # Representation choice
//!
//! In this module, the curve parameters are stored as the pair $(A, B)$.
//! The native point representation for arithmetic is **x-only projective**
//! coordinates on the Kummer line, rather than full affine points.
//!
//! # Why x-only arithmetic?
//!
//! On a Montgomery curve, scalar multiplication can be implemented using only
//! $x$-coordinates via the Montgomery ladder. This works on the Kummer quotient
//! $E / \{\pm 1\}$, where a point $P$ and its inverse $-P$ have the same image.
//!
//! The advantage is that the ladder uses a uniform sequence of differential
//! additions and doublings, which is especially convenient for constant-time
//! scalar multiplication.

<<<<<<< HEAD
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use fp::field_ops::FieldOps;
use fp::{ref_field_impl, ref_field_trait_impl};
=======
use core::fmt;
use fp::field_ops::{FieldOps, FieldRandom};

>>>>>>> origin/main
use crate::curve_ops::Curve;
use crate::point_montgomery::KummerPoint;

/// A Montgomery curve
///
/// $By^2 = x(x^2+ Ax + 1)$
///
/// over a field `F`.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct MontgomeryCurve<F: FieldOps> {
    /// The invariant $A$ in the equation of the montgomery curve
    pub a: F,
    /// The invariant $B$ in the equation of the montgomery curve
    pub b: F,
}

<<<<<<< HEAD
ref_field_impl! {
    impl<F> MontgomeryCurve<F> {

        pub fn new(a: F, b: F) -> Self {
            assert!(Self::is_smooth(&a, &b));
            Self { a, b }
=======

impl<F> fmt::Display for MontgomeryCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if F::characteristic()[0] != 2 {
            if f.alternate() {
                write!(
                    f,
                    "MontgomeryCurve {{\n  By^2 = x(x^2 + Ax + 1)\n  A = {}\n  B = {}\n}}",
                    self.a, self.b
                )
            } else {
                write!(
                    f,
                    "By^2 = x(x^2 + Ax + 1) over (A={}, B={})",
                    self.a, self.b
                )
            }
        } else {
            if f.alternate() {
                write!(
                    f,
                    "MontgomeryCurve {{\n  y^2 + xy = x(x^2 + Ax + B^2)\n  A = {}\n  B = {}\n}}",
                    self.a, self.b
                )
            } else {
                write!(
                    f,
                    "y^2 + xy = x(x^2 + Ax + B^2) over (A={}, B={})",
                    self.a, self.b
                )
            }
        }
    }
}
impl<F: FieldOps + Copy> MontgomeryCurve<F> {
    // -------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------

    /// Construct a Montgomery curve from its two coefficients `A` and `B`.
    pub fn new(a: F, b: F) -> Self {
        assert!(Self::is_smooth(&a, &b));
        Self { a, b }
    }

    /// True if and only if the curve is smooth
    pub fn is_smooth(a: &F, b: &F) -> bool {
        if F::characteristic()[0] != 2 {
            // In odd char, smoothness criterion is just b != 0 and a != 2, -2
            let two = F::one().double();
            !bool::from(b.is_zero()) && *a != two && *a != two.negate()
        } else {
            // In char 2, smoothness criterion is just b != 0
            !bool::from(b.is_zero())
>>>>>>> origin/main
        }

<<<<<<< HEAD
        pub fn is_smooth(a: &F, b: &F) -> bool {
            if F::characteristic()[0] != 2 {
                let two = <F as FieldOps>::double(&F::one());
                let minus_two = -&two;

                !bool::from(b.is_zero()) && a != &two && a != &minus_two
            } else {
                !bool::from(b.is_zero())
            }
        }
=======
    // -------------------------------------------------------------------
    // Invariants and constants
    // -------------------------------------------------------------------
>>>>>>> origin/main

        pub fn a_invariants(&self) -> [F; 2] {
            [self.a.clone(), self.b.clone()]
        }

        pub fn a24(&self) -> F {
            assert_ne!(F::characteristic()[0], 2);

            let two = <F as FieldOps>::double(&F::one());
            let four = <F as FieldOps>::double(&two);
            let four_inv = <F as FieldOps>::invert(&four).unwrap();

            let tmp = &self.a + &two;
            &tmp * &four_inv
        }
    }
}

impl<F: FieldOps + Copy + FieldRandom> MontgomeryCurve<F> {
    /// Sample a random Kummer/x-line point on this curve using the provided RNG.
    ///
    /// The method repeatedly samples an affine `x`-coordinate and returns the
    /// corresponding projective Kummer point `(x:1)` once it represents a valid
    /// point on the curve.
    pub fn random_point(
        &self,
        rng: &mut (impl rand::CryptoRng + rand::Rng),
    ) -> KummerPoint<F> {
        loop {
            let x = F::random(rng);
            let p = KummerPoint::from_x(x);
            if self.is_on_curve(&p) {
                return p;
            }
        }
    }
}


// -------------------------------------------------------------------
// Curve predicates
// -------------------------------------------------------------------

<<<<<<< HEAD
ref_field_trait_impl! {
    impl<F> Curve for MontgomeryCurve<F> {

        type BaseField = F;
        type Point = KummerPoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            if point.is_identity() {
                return true;
            }

            let z_inv = <F as FieldOps>::invert(&point.z);
            if bool::from(z_inv.is_none()) {
                return true;
            }

            let x = &point.x * &z_inv.unwrap();

            if F::characteristic()[0] != 2 {
                let xsq = <F as FieldOps>::square(&x);
                let xcubed = &x * &xsq;
                let axsq = &self.a * &xsq;

                let sum1 = &xcubed + &axsq;
                let sum = &sum1 + &x;

                let binv = <F as FieldOps>::invert(&self.b).unwrap();
                let val = &sum * &binv;

                <F as FieldOps>::legendre(&val) > -1i8
=======
impl<F: FieldOps + Copy+ FieldRandom> Curve for MontgomeryCurve<F> {
    type BaseField = F;
    type Point = KummerPoint<F>;

    /// Return `true` if `point` is a valid Kummer/x-line point for this curve.
    fn is_on_curve(&self, point: &Self::Point) -> bool {
        if point.is_identity() {
            true
        } else {
            // Normalise projective (X : Z) to affine x = X / Z.
            let z_inv = <F as FieldOps>::invert(&point.z);
            if bool::from(z_inv.is_none()) {
                // Z = 0 means identity, already handled above.
                return true;
            }
            let x = point.x * z_inv.unwrap();
            if F::characteristic()[0] != 2 {
                // Check if `(x³ + A x² + x)/B` is a square in F or not.
                let xsq = <F as FieldOps>::square(&x);
                let xcubed = x * xsq;
                let axsq = self.a * xsq;
                let binv = <F as FieldOps>::invert(&self.b).unwrap();
                let sum = xcubed + axsq + x;
                <F as FieldOps>::legendre(&(sum * binv)) > -1i8
>>>>>>> origin/main
            } else {
                if bool::from(x.is_zero()) {
                    true
                } else {
<<<<<<< HEAD
=======
                    // If x != 0, set t = y/x, which gives the equation t^2 + t = x + a + b^2/x.
                    // So x is the x-coordinate of a point on the curve iff there is some t in F
                    // s.t. t^2 + t = x + a + b^2/x, which translates into
                    // Tr_{F/F_2}(x + a + b^2/x) == 0.
>>>>>>> origin/main
                    let bsq = <F as FieldOps>::square(&self.b);
                    let xinv = x.invert().unwrap();

                    let tmp1 = &x + &self.a;
                    let tmp2 = &bsq * &xinv;
                    let rhs = &tmp1 + &tmp2;

                    bool::from(rhs.trace().is_zero())
                }
            }
        }

<<<<<<< HEAD
        fn random_point(&self) -> Self::Point {
            todo!()
        }

        fn j_invariant(&self) -> F {
            if F::characteristic()[0] != 2 {
                let two = <F as FieldOps>::double(&F::one());
                let three = &two + &F::one();
                let four = <F as FieldOps>::double(&two);
                let sixteen = <F as FieldOps>::square(&four);
                let twofivesix = <F as FieldOps>::square(&sixteen);

                let asq = <F as FieldOps>::square(&self.a);
                let asq_min_three = &asq - &three;

                let asq_min_three_sq = <F as FieldOps>::square(&asq_min_three);
                let asq_min_three_cubed = &asq_min_three * &asq_min_three_sq;

                let denom = &asq - &four;
                let denom_inv = <F as FieldOps>::invert(&denom)
                    .into_option()
                    .expect("a != ±2");

                let tmp = &twofivesix * &asq_min_three_cubed;
                &tmp * &denom_inv
            } else {
                let bsq = <F as FieldOps>::square(&self.b);
                let bfourth = <F as FieldOps>::square(&bsq);
                <F as FieldOps>::invert(&bfourth).unwrap()
            }
        }

        fn a_invariants(&self) -> Vec<Self::BaseField> {
            MontgomeryCurve::a_invariants(self).to_vec()
        }
    }
}



// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F> ConditionallySelectable for MontgomeryCurve<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self{
            a: F::conditional_select(& a.a, &b.a, choice),
            b: F::conditional_select(& a.b, &b.b, choice),
=======
    fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
        MontgomeryCurve::random_point(self, rng)
    }

    fn j_invariant(&self) -> F {
        if F::characteristic()[0] != 2 {
            // The j-invariant of `B y² = x(x² + A x + 1)` is 256*(a² - 3)³/(a² - 4)
            let two = <F as FieldOps>::double(&F::one());
            let three = two + F::one();
            let four = <F as FieldOps>::double(&two);
            let sixteen = <F as FieldOps>::square(&four);
            let twofivesix = <F as FieldOps>::square(&sixteen);

            let asq = <F as FieldOps>::square(&self.a);
            let asq_min_three = asq - three;
            let asq_min_three_cubed = asq_min_three * <F as FieldOps>::square(&asq_min_three);

            let asq_min_four_inv = <F as FieldOps>::invert(&(asq - four))
                .into_option()
                .expect("a should be different from 2, -2.");

            twofivesix * asq_min_three_cubed * asq_min_four_inv
        } else {
            // The j-invariant of `y² + xy = x(x² + A x + B²)` is 1/b⁴
            assert!(!bool::from(self.b.is_zero()));
            let bsq = <F as FieldOps>::square(&self.b);
            let bfourth = <F as FieldOps>::square(&bsq);
            <F as FieldOps>::invert(&bfourth).unwrap()
>>>>>>> origin/main
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        F::conditional_assign(&mut self.a, &other.a, choice);
        F::conditional_assign(&mut self.b, &other.b, choice);
    }
<<<<<<< HEAD

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.a, &mut b.a, choice);
        F::conditional_swap(&mut a.b, &mut b.b, choice);
    }
}

impl<F> ConstantTimeEq for MontgomeryCurve<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.a.ct_eq(&other.a) & self.b.ct_eq(&other.b)
    }

    fn ct_ne(&self, other: &Self) -> Choice { !self.ct_eq(other) }
}



=======
}
>>>>>>> origin/main
