//! Elliptic curve definition in Edwards form.
//!
//! # Odd characteristic
//!
//! The curve is given by
//!
//! $$
//! x^2 + y^2 = 1 + d x^2 y^2
//! $$
//!
//! Parameters: $d_1 = 0$ (unused), $d_2 = d$.
//! Identity: $(0, 1)$. Complete when $d$ is not a square.
//!
//! # Characteristic $2$
//!
//! The curve is given by
//!
//! $$
//! d_1(x + y) + d_2(x^2 + y^2)
//! = xy + xy(x + y) + x^2 y^2
//! $$
//!
//! Parameters: $d_1, d_2$ with $d_1 \ne 0$, $d_2 \ne d_1^2 + d_1$.
//! Identity: $(0, 0)$. Complete when $\mathrm{Tr}(d_2) = 1$.
//!
//! # References
//!
//! - Binary Edwards curves (characteristic $2$): Bernstein–Lange–Rezaeian Farashahi (2008)
//! - Odd characteristic: <https://hyperelliptic.org/EFD/g1p/auto-edwards.html>

use core::fmt;
use fp::field_ops::{FieldOps, FieldRandom};

use crate::curve_ops::Curve;
use crate::point_edwards::EdwardsPoint;

/// An Edwards curve over a field `F`, covering both odd and even characteristic.
///
/// In odd characteristic only `d2` is used (the parameter `d`).
/// In characteristic 2 both `d1` and `d2` are used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdwardsCurve<F: FieldOps> {
    /// The invariant d1 in the equation
    pub d1: F,
    /// The invariant d2 in the equation
    pub d2: F,
}

impl<F> fmt::Display for EdwardsCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if F::characteristic()[0] != 2 {
            if f.alternate() {
                write!(
                    f,
                    "EdwardsCurve {{\n  x^2 + y^2 = 1 + d x^2 y^2\n  d = {}\n}}",
                    self.d2
                )
            } else {
                write!(f, "x^2 + y^2 = 1 + ({})x^2y^2", self.d2)
            }
        } else {
            if f.alternate() {
                write!(
                    f,
                    "EdwardsCurve {{\n  d1(x+y) + d2(x^2+y^2) = xy + xy(x+y) + x^2y^2\n  d1 = {}\n  d2 = {}\n}}",
                    self.d1, self.d2
                )
            } else {
                write!(
                    f,
                    "{}(x+y) + {}(x^2+y^2) = xy + xy(x+y) + x^2y^2",
                    self.d1, self.d2
                )
            }
        }
    }
}

impl<F: FieldOps + FieldRandom> EdwardsCurve<F> {
    /// Construct an odd-characteristic Edwards curve `x² + y² = 1 + d·x²·y²`.
    ///
    /// Stores `d` as `d2`; `d1` is set to zero (unused).
    pub fn new(d: F) -> Self {
        assert!(F::characteristic()[0] != 2, "use new_binary() for char 2");
        assert!(d != F::zero(), "d must be nonzero");
        assert!(d != F::one(), "d must not be 1");
        Self {
            d1: F::zero(),
            d2: d,
        }
    }

    /// Construct a binary Edwards curve
    /// `d₁(x+y) + d₂(x²+y²) = xy + xy(x+y) + x²y²`.
    pub fn new_binary(d1: F, d2: F) -> Self {
        assert_eq!(F::characteristic()[0], 2, "use new() for odd char");
        assert!(d1 != F::zero(), "d1 must be nonzero");
        let d1_sq = <F as FieldOps>::square(&d1);
        assert!(d2 != d1_sq + d1, "d2 must differ from d1² + d1");
        Self { d1, d2 }
    }

    /// Convenience accessor: the parameter `d` in odd characteristic.
    pub fn d(&self) -> F {
        self.d2
    }

    /// Check whether the affine point `(x, y)` lies on the curve.
    pub fn contains(&self, x: &F, y: &F) -> bool {
        if F::characteristic()[0] != 2 {
            // x² + y² == 1 + d·x²·y²
            let x2 = <F as FieldOps>::square(x);
            let y2 = <F as FieldOps>::square(y);
            x2 + y2 == F::one() + self.d2 * x2 * y2
        } else {
            // d₁(x+y) + d₂(x²+y²) == xy + xy(x+y) + x²y²
            let x2 = <F as FieldOps>::square(x);
            let y2 = <F as FieldOps>::square(y);
            let xy = *x * *y;
            let xpy = *x + *y;
            self.d1 * xpy + self.d2 * (x2 + y2) == xy + xy * xpy + x2 * y2
        }
    }

    /// Sample a random affine point on this Edwards curve using the provided RNG.
    ///
    /// This currently uses a square-root-based construction and is implemented
    /// only for odd characteristic. It returns a point `P` such that
    /// `self.is_on_curve(&P)` holds.
    pub fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> EdwardsPoint<F> {
        assert!(
            F::characteristic()[0] != 2,
            "random_point currently implemented only for odd characteristic"
        );

        loop {
            let x = F::random(rng);
            let x2 = <F as FieldOps>::square(&x);
            let denom = F::one() - self.d2 * x2;

            if bool::from(denom.is_zero()) {
                continue;
            }

            let rhs = (F::one() - x2) * denom.invert().unwrap();

            if let Some(y) = rhs.sqrt().into_option() {
                let p = EdwardsPoint::new(x, y);
                debug_assert!(self.is_on_curve(&p));
                return p;
            }
        }
    }
}

impl<F: FieldOps + FieldRandom> Curve for EdwardsCurve<F> {
    type BaseField = F;
    type Point = EdwardsPoint<F>;

    fn is_on_curve(&self, point: &Self::Point) -> bool {
        self.contains(&point.x, &point.y)
    }

    fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
       EdwardsCurve::random_point(self, rng)
    }

    fn j_invariant(&self) -> F {
        if F::characteristic()[0] != 2 {
            // j = 16(1 + 14d + d²)³ / (d(1 − d)⁴)
            let d = self.d2;
            let d2 = <F as FieldOps>::square(&d);
            let two = <F as FieldOps>::double(&F::one());
            let four = <F as FieldOps>::double(&two);
            let eight = <F as FieldOps>::double(&four);
            let fourteen = eight + four + two;
            let sixteen = <F as FieldOps>::double(&eight);

            let inner = F::one() + fourteen * d + d2;
            let inner_cubed = inner * <F as FieldOps>::square(&inner);
            let numer = sixteen * inner_cubed;

            let one_minus_d = F::one() - d;
            let omd2 = <F as FieldOps>::square(&one_minus_d);
            let omd4 = <F as FieldOps>::square(&omd2);
            let denom = d * omd4;

            numer
                * denom
                    .invert()
                    .into_option()
                    .expect("d(1-d)^4 must be invertible")
        } else {
            // j = 1 / (d₁⁴ (d₁⁴ + d₁² + d₂²))
            let d1_sq = <F as FieldOps>::square(&self.d1);
            let d1_4 = <F as FieldOps>::square(&d1_sq);
            let d2_sq = <F as FieldOps>::square(&self.d2);
            let denom = d1_4 * (d1_4 + d1_sq + d2_sq);
            denom
                .invert()
                .into_option()
                .expect("j-invariant denominator must be invertible")
        }
    }

    fn a_invariants(&self) -> Vec<Self::BaseField> {
        if F::characteristic()[0] != 2 {
            vec![self.d2]
        } else {
            vec![self.d1, self.d2]
        }
    }
}
