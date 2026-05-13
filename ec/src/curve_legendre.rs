//! Elliptic curve definition in Legendre form.
//!
//! # Equation
//!
//! The **Legendre** family over a field `F` of odd characteristic is
//!
//! $$
//! E_\lambda : y^2 = x(x-1)(x-\lambda),
//! \qquad \lambda \neq 0,1.
//! $$
//!
//! Expanding the right-hand side gives the equivalent Weierstrass model
//!
//! $$
//! y^2 = x^3 - (1+\lambda)x^2 + \lambda x,
//! $$
//!
//! so the associated $a$-invariants are
//!
//! $$
//! (a_1,a_2,a_3,a_4,a_6) = (0,\,-(1+\lambda),\,0,\,\lambda,\,0).
//! $$
//!
//! # Distinguished points
//!
//! The Legendre form makes the full rational $2$-torsion visible:
//!
//! $$
//! O,\quad (0,0),\quad (1,0),\quad (\lambda,0).
//! $$
//!
//! The group identity is the point at infinity `O`.
//!
//! # j-invariant
//!
//! The $j$-invariant of `E_λ` is
//!
//! $$
//! j(E_\lambda)=256\frac{(\lambda^2-\lambda+1)^3}{\lambda^2(\lambda-1)^2}.
//! $$
//!
//! # References
//!
//! - Roland Auer and Jaap Top, *Legendre Elliptic Curves over Finite Fields*,
//!   Journal of Number Theory 95 (2002), 303–312.
//! - HongFeng Wu and RongQuan Feng,
//!   *On the isomorphism classes of Legendre elliptic curves over finite fields*,
//!   Sci. China Math. 54(9) (2011), 1885–1890.
use crate::curve_ops::Curve;
use crate::curve_weierstrass::WeierstrassCurve;
use crate::point_legendre::LegendrePoint;
use core::fmt;
use fp::field_ops::{FieldOps, FieldRandom};
use fp::{ref_field_impl, ref_field_trait_impl};

/// A Legendre elliptic curve over a field `F`.
///
/// This stores the parameter `λ` of the curve
///
/// $$
/// E_\lambda : y^2 = x(x-1)(x-\lambda).
/// $$
///
/// The nonsingular case is exactly `λ != 0` and `λ != 1`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegendreCurve<F: FieldOps> {
    /// The Legendre parameter $\lambda$ in:
    /// $$
    /// E_{\lambda} : y^2 = x(x - 1)(x - \lambda)
    /// $$
    ///
    /// For a nonsingular Legendre curve one must have $\lambda \neq 0$ and $\lambda \neq 1$.
    pub lambda: F,
}

ref_field_impl! {
    impl<F: FieldOps + FieldRandom> LegendreCurve<F> {
        /// Construct the Legendre curve $y^2 = x(x-1)(x-\lambda).$
        pub fn new(lambda: F) -> Self {
            Self { lambda }
        }

        /// Returns `true` if and only if the model is singular.
        ///
        /// For Legendre form, singularity occurs exactly when `λ ∈ {0,1}`
        pub fn is_singular(&self) -> bool {
            self.lambda == F::zero() || self.lambda == F::one()
        }

        /// Returns the affine right-hand side
        ///
        /// $$
        /// x(x-1)(x-\lambda).
        /// $$
        pub fn rhs(&self, x: &F) -> F {
            x * &(&(x - &F::one()) * &(x - &self.lambda))
        }

        /// Checks whether the affine point `(x, y)` satisfies
        ///
        /// $$
        /// y^2 = x(x-1)(x-\lambda).
        /// $$
        pub fn contains(&self, x: &F, y: &F) -> bool {
            let lhs = <F as FieldOps>::square(y);
            let rhs = self.rhs(x);
            lhs == rhs
        }

        /// Samples a random affine point on the curve.
        ///
        /// This uses a square-root based strategy in odd characteristic:
        /// choose a random `x`, compute `rhs(x)`, and return `(x, y)` whenever
        /// `rhs(x)` is a square.
        ///
        /// # Panics
        ///
        /// Panics if the curve is singular.
        pub fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> LegendrePoint<F> {
            assert!(
                !self.is_singular(),
                "cannot sample points on a singular Legendre curve"
            );

            loop {
                let x = F::random(rng);
                let rhs = self.rhs(&x);

                if let Some(y) = rhs.sqrt().into_option() {
                    let use_neg = (rng.next_u32() & 1) == 1;
                    let y = if use_neg { -&y } else { y };
                    return LegendrePoint::new(x, y);
                }
            }
        }

        /// Returns the Weierstrass `j`-invariant of `E_λ`
        ///
        /// j(E_\lambda)=256\frac{(\lambda^2-\lambda+1)^3}{\lambda^2(\lambda-1)^2}.
        ///
        pub fn j_invariant_model(&self) -> F {
            assert!(
                !self.is_singular(),
                "j-invariant is undefined for a singular Legendre curve"
            );

            let lambda_sq = <F as FieldOps>::square(&self.lambda);
            let t = &(&lambda_sq - &self.lambda) + &F::one();
            let tsq = <F as FieldOps>::square(&t);
            let num = &F::from_u64(256) * &(&t * &tsq);

            let lambda_minus_one = &self.lambda - &F::one();
            let den = &lambda_sq * &<F as FieldOps>::square(&lambda_minus_one);

            let den_inv = den
                .invert()
                .into_option()
                .expect("denominator must be invertible for λ != 0,1");

            &num * &den_inv
        }

        /// Returns the a-invariants of
        ///
        /// y² = x³ - (1+λ)x² + λx
        ///
        /// i.e. [a1, a2, a3, a4, a6] = [0, -(1+λ), 0, λ, 0].
        pub fn a_invariants(&self) -> [F; 5] {
            [
                F::zero(),
                -&(&F::one() + &self.lambda),
                F::zero(),
                self.lambda,
                F::zero(),
            ]
        }

        /// Returns the same curve written in general Weierstrass form.
        ///
        /// The Legendre equation
        ///
        /// `y² = x(x-1)(x-λ)`
        ///
        /// expands to
        ///
        /// `y² = x³ - (1+λ)x² + λx`,
        ///
        /// so the corresponding Weierstrass coefficients are
        ///
        /// `[a1, a2, a3, a4, a6] = [0, -(1+λ), 0, λ, 0]`.
        pub fn to_weierstrass(&self) -> WeierstrassCurve<F> {
            WeierstrassCurve::new(
                F::zero(),
                -&(&F::one() + &self.lambda),
                F::zero(),
                self.lambda,
                F::zero(),
            )
        }

        /// Returns a short Weierstrass model isomorphic to this Legendre curve.
        ///
        /// Over fields of characteristic different from `2` and `3`, the change of
        /// variables
        ///
        /// `x = X + (1+λ)/3`, `y = Y`
        ///
        /// transforms
        ///
        /// `y² = x³ - (1+λ)x² + λx`
        ///
        /// into
        ///
        /// `Y² = X³ + AX + B`
        ///
        /// where
        ///
        /// `A = -(λ² - λ + 1)/3`
        ///
        /// `B = -((λ + 1)(λ - 2)(2λ - 1))/27`.
        ///
        /// # Panics
        ///
        /// Panics if the characteristic is `<= 3`, or if the curve is singular.
        pub fn to_short_weierstrass(&self) -> WeierstrassCurve<F> {
            assert!(
                !self.is_singular(),
                "cannot convert a singular Legendre curve to short Weierstrass form"
            );
            assert!(
                F::characteristic()[0] > 3,
                "short Weierstrass conversion requires characteristic != 2, 3"
            );

            let three = F::from_u64(3);
            let twenty_seven = F::from_u64(27);

            let three_inv = three
                .invert()
                .into_option()
                .expect("3 must be invertible in characteristic != 3");

            let twenty_seven_inv = twenty_seven
                .invert()
                .into_option()
                .expect("27 must be invertible in characteristic != 3");

            let lambda_sq = <F as FieldOps>::square(&self.lambda);
            let lmsq_minus_lm = &lambda_sq - &self.lambda;

            let num_a = -&(&lmsq_minus_lm + &F::one());
            let a = &num_a * &three_inv;

            let lm_plus_one = &self.lambda + &F::one();
            let lm_minus_two = &self.lambda - &F::from_u64(2);
            let twolm_minus_one = &<F as FieldOps>::double(&self.lambda) - &F::one();

            let num_b = -&(&lm_plus_one * &(&lm_minus_two * &twolm_minus_one));
            let b = &num_b * &twenty_seven_inv;

            WeierstrassCurve::new_short(a, b)
        }
    }
}

impl<F> fmt::Display for LegendreCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(
                f,
                "LegendreCurve {{\n  y^2 = x(x-1)(x - {})\n}}",
                self.lambda
            )
        } else {
            write!(f, "y^2 = x(x-1)(x - {})", self.lambda)
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + FieldRandom> Curve for LegendreCurve<F> {
        type BaseField = F;
        type Point = LegendrePoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            if point.infinity {
                true
            } else {
                let lhs = <F as FieldOps>::square(&point.y);
                let rhs = self.rhs(&point.x);
                lhs == rhs
            }
        }

        fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
            LegendreCurve::random_point(self, rng)
        }

        /// Returns the $j$-invariant of the curve.
        ///
        /// This is a complete invariant of elliptic curves over algebraically
        /// closed fields up to isomorphism.
        fn j_invariant(&self) -> F {
            LegendreCurve::j_invariant_model(&self)
        }

        /// Returns the $a$-invariants as a vector.
        fn a_invariants(&self) -> Vec<Self::BaseField> {
            LegendreCurve::a_invariants(self).to_vec()
        }
    }
}
