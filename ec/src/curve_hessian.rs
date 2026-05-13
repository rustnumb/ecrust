//! Elliptic curve definition in (generalized) Hessian form.
//!
//! # Equation
//!
//! We use the generalized Hessian model from Farashahi--Joye:
//!
//! $$
//! x^3 + y^3 + c = dxy,
//! $$
//!
//! together with its projective closure
//!
//! $$
//! X^3 + Y^3 + c Z^3 = d X Y Z.
//! $$
//!
//! The neutral element is the point at infinity
//! $$O = (1 : -1 : 0).$$
//!
//! ## Parameter conventions
//!
//! The paper `Efficient Arithmetic on Hessian Curves` uses the parameterization
//! $x^3 + y^3 + c = dxy$, while the EFD page for ordinary Hessian curves uses
//! $x^3 + y^3 + 1 = 3 d_{e} xy$.
//!
//! In this module we store the paper parameter $d$, so the relation is
//!
//! $d_{p} = 3 * d_{e}$.
//!
//! For ordinary Hessian curves, set $c = 1$.
//!
//! # Smoothness
//!
//! The generalized Hessian curve is nonsingular when
//!
//! $$c \neq 0 $$
//! $$d^3 \neq 27c.$$
//!
//! # References
//!
//! - Reza R. Farashahi and Marc Joye,
//!   *Efficient Arithmetic on Hessian Curves*, PKC 2010.
//! - Explicit-Formulas Database (EFD), Hessian curves:
//!   <https://www.hyperelliptic.org/EFD/g1p/auto-hessian.html>
//! - EFD, projective Hessian formulas:
//!   <https://www.hyperelliptic.org/EFD/g1p/auto-hessian-standard.html>

use crate::curve_ops::Curve;
use crate::curve_weierstrass::WeierstrassCurve;
use crate::point_hessian::HessianPoint;
use crate::point_weierstrass::AffinePoint;
use core::fmt;
use fp::field_ops::{FieldOps, FieldRandom};
use fp::{ref_field_impl, ref_field_trait_impl};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// A generalized Hessian curve
///
/// $$x^3 + y^3 + c = dxy.$$
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct HessianCurve<F: FieldOps> {
    /// The constant term `c` in `x^3 + y^3 + c = dxy`.
    pub c: F,
    /// The paper parameter `d` in `x^3 + y^3 + c = dxy`.
    pub d: F,
}

impl<F> fmt::Display for HessianCurve<F>
where
    F: FieldOps + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(
                f,
                "HessianCurve {{\n  x^3 + y^3 + c = dxy\n  c = {}\n  d = {}\n}}",
                self.c, self.d
            )
        } else {
            write!(f, "x^3 + y^3 + ({}) = ({})xy", self.c, self.d)
        }
    }
}

ref_field_impl! {
    impl<F: FieldOps + FieldRandom> HessianCurve<F> {
        /// Construct a generalized Hessian curve
        ///
        /// $$x^3 + y^3 + c = dxy.$$
        pub fn new(c: F, d: F) -> Self {
            assert!(Self::is_smooth(&c, &d), "singular Hessian curve");
            Self { c, d }
        }

        /// Construct an ordinary Hessian curve
        ///
        /// $$x^3 + y^3 + 1 = dxy.$$
        pub fn new_hessian(d: F) -> Self {
            Self::new(F::one(), d)
        }

        /// Construct an ordinary Hessian curve from the EFD parameterization
        ///
        /// $$x^3 + y^3 + 1 = 3d_{efd}xy.$$
        ///
        /// Internally we store `d = 3*d_efd` to match the paper's convention.
        pub fn new_efd(d_efd: F) -> Self {
            Self::new(F::one(),  &F::from_u64(3) * &d_efd)
        }

        /// Return `true` if the generalized Hessian discriminant is nonzero.
        ///
        /// The smoothness criterion is
        ///
        /// $$c \neq 0\;\text{and}\; d^3 \neq 27c.$$
        pub fn is_smooth(c: &F, d: &F) -> bool {
            if bool::from(c.is_zero()) {
                return false;
            }

            let d2 = <F as FieldOps>::square(d);
            let d3 = d * &d2;
            d3 !=  &F::from_u64(27) * c
        }

        /// Check whether the affine point `(x, y)` satisfies
        ///
        /// $$x^3 + y^3 + c = dxy.$$
        pub fn contains_affine(&self, x: &F, y: &F) -> bool {
            let x2 = <F as FieldOps>::square(x);
            let y2 = <F as FieldOps>::square(y);
            let x3 = x * &x2;
            let y3 = y * &y2;
            &x3 + &(&y3 + &self.c) == &self.d * &(x * y)
        }

        /// Check whether the projective point `(X:Y:Z)` satisfies
        ///
        /// $$X^3 + Y^3 + c Z^3 = dXYZ.$$
        pub fn contains_projective(&self, x: &F, y: &F, z: &F) -> bool {
            let x2 = <F as FieldOps>::square(x);
            let y2 = <F as FieldOps>::square(y);
            let z2 = <F as FieldOps>::square(z);
            let x3 = x * &x2;
            let y3 = y * &y2;
            let z3 = z * &z2;
            &x3 + &(&y3 + &(&self.c * &z3)) == &self.d * &(x * &(y * z))
        }

        /// Return `[c, d]`.
        pub fn a_invariants(&self) -> [F; 2] {
            [self.c, self.d]
        }

        /// Return the neutral element `(1:-1:0)`.
        pub fn neutral_point(&self) -> HessianPoint<F> {
            HessianPoint::identity()
        }

        /// Best-effort random point sampling.
        ///
        /// A generic implementation of cube solving is not yet available in the
        /// field layer, so this method performs a bounded affine rejection search
        /// and falls back to the neutral element if no finite point is found.
        pub fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> HessianPoint<F> {
            loop {
                let x = F::random(rng);
                let y = F::random(rng);
                if self.contains_affine(&x, &y) {
                    let p = HessianPoint::from_affine(x, y);
                    debug_assert!(self.is_on_curve(&p));
                    return p;
                }
            }
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + FieldRandom> Curve for HessianCurve<F> {
        type BaseField = F;
        type Point = HessianPoint<F>;

        fn is_on_curve(&self, point: &Self::Point) -> bool {
            self.contains_projective(&point.x, &point.y, &point.z)
        }

        fn random_point(&self, rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self::Point {
            HessianCurve::random_point(self, rng)
        }

        fn j_invariant(&self) -> F {
            // Farashahi--Joye, Eq. (2) / Appendix A:
            // j(H_{c,d}) = c^{-1} * [ d(d^3 + 216c)/(d^3 - 27c) ]^3.
            let d2 = <F as FieldOps>::square(&self.d);
            let d3 = &self.d * &d2;
            let twenty_seven_c =  &F::from_u64(27) * &self.c;
            let two_hundred_sixteen_c = &F::from_u64(216) * &self.c;

            let inner_num = &self.d * &(&d3 + &two_hundred_sixteen_c);
            let inner_den = &d3 - &twenty_seven_c;
            let inner = &inner_num
                * &inner_den
                    .invert()
                    .into_option()
                    .expect("Hessian j-invariant denominator must be invertible");
            let inner_cubed = &inner * &<F as FieldOps>::square(&inner);

            &self.c
                .invert()
                .into_option()
                .expect("c must be invertible on a nonsingular Hessian curve")
                * &inner_cubed
        }

        fn a_invariants(&self) -> Vec<Self::BaseField> {
            HessianCurve::a_invariants(self).to_vec()
        }
    }
}

impl<F> ConditionallySelectable for HessianCurve<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            c: F::conditional_select(&a.c, &b.c, choice),
            d: F::conditional_select(&a.d, &b.d, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.c.conditional_assign(&other.c, choice);
        self.d.conditional_assign(&other.d, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        F::conditional_swap(&mut a.c, &mut b.c, choice);
        F::conditional_swap(&mut a.d, &mut b.d, choice);
    }
}

impl<F> ConstantTimeEq for HessianCurve<F>
where
    F: FieldOps + Copy + ConstantTimeEq,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c.ct_eq(&other.c) & self.d.ct_eq(&other.d)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}

ref_field_impl! {
    impl<F: FieldOps> HessianCurve<F> {
        /// Return the EFD Hessian parameter `delta` for the ordinary Hessian model
        ///
        ///   $$x'^3 + y'^3 + 1 = 3*delta*x'*y'$$
        ///
        /// obtained from the generalized Hessian
        ///
        ///   $x^3 + y^3 + c = d*x*y$
        ///
        /// by $x' = x/zeta$, $y' = y/zeta$ with $zeta^3 = c$.
        ///
        /// Since this crate stores the paper parameter `d`, the relation is
        ///
        ///   $delta = d / (3*zeta)$.
        pub fn efd_parameter_with_zeta(&self, zeta: F) -> Option<F> {
            let zeta_sq = <F as FieldOps>::square(&zeta);
            let zeta_cu = &zeta * &zeta_sq;
            if zeta_cu != self.c {
                return None;
            }

            let denom =  &F::from_u64(3) * &zeta;
            denom.invert().into_option().map(|inv| &self.d * &inv)
        }

        /// Convert this generalized Hessian curve to the EFD short Weierstrass model.
        ///
        /// EFD for $x^3 + y^3 + 1 = 3*delta*x*y$ gives:
        ///
        ///   $$v^2 = u^3 + a*u + b$$
        ///
        /// with
        ///
        ///   $$a = -27*delta*(delta^3 + 8)$$
        ///   $$b =  54*(delta^6 - 20*delta^3 - 8)$$
        ///
        /// This returns `None` if `zeta^3 != c` or if the map is unavailable
        /// in the current characteristic (e.g. 2 or 3).
        pub fn to_weierstrass_curve_with_zeta(&self, zeta: F) -> Option<WeierstrassCurve<F>> {
            let delta = self.efd_parameter_with_zeta(zeta)?;

            let delta_sq = <F as FieldOps>::square(&delta);
            let delta_cu = &delta * &delta_sq;
            let delta_six = &delta_cu * &delta_cu;

            let delta_cu_pl_8 = &delta_cu + &F::from_u64(8);
            let tmp1 = &delta * &delta_cu_pl_8;
            let a = -&(&F::from_u64(27) * &tmp1);

            let twenty_delta_cu = &F::from_u64(20) * &delta_cu;
            let sum1 = &delta_six - &twenty_delta_cu;
            let sum2 = &sum1 - &F::from_u64(8);

            let b = &F::from_u64(54) * &sum2;

            Some(WeierstrassCurve {
                a1: F::zero(),
                a2: F::zero(),
                a3: F::zero(),
                a4: a,
                a6: b,
            })
        }

        /// Map a Hessian point to the EFD short Weierstrass model.
        ///
        /// For the ordinary Hessian x'^3 + y'^3 + 1 = 3*delta*x'*y',
        /// the EFD birational map is
        ///
        ///   $$u = 12*(delta^3 - 1)/(delta + x' + y') - 9*delta^2$$
        ///   $$v = 36*(y' - x')*(delta^3 - 1)/(delta + x' + y')$$
        ///
        /// For the generalized Hessian we first rescale by zeta:
        ///
        ///   $$x' = x/zeta$$
        /// , $$y' = y/zeta$$,
        /// $$zeta^3 = c$$.
        pub fn map_point_to_weierstrass_with_zeta(
            &self,
            p: &HessianPoint<F>,
            zeta: F,
        ) -> Option<AffinePoint<F>> {
            if p.is_identity() {
                return Some(AffinePoint::identity());
            }

            let (x, y) = p.to_affine()?;
            let zeta_inv = zeta.invert().into_option()?;
            let xh = &x * &zeta_inv;
            let yh = &y * &zeta_inv;

            let delta = self.efd_parameter_with_zeta(zeta)?;
            let delta_sq = <F as FieldOps>::square(&delta);
            let delta_cu = &delta * &delta_sq;
            let dm1 = &delta_cu - &F::one();

            let denom = &delta + &(&xh + &yh);
            let denom_inv = match denom.invert().into_option() {
                Some(v) => v,
                None => return Some(AffinePoint::identity()),
            };

            let u =  &(&F::from_u64(12) * &(&dm1 * &denom_inv)) - &(&F::from_u64(9) * &delta_sq);
            let v =  &F::from_u64(36) * &(&(&yh - &xh) * &(&dm1 * &denom_inv));

            Some(AffinePoint::new(u, v))
        }

        /// Inverse EFD birational map from the short Weierstrass model back to
        /// the generalized Hessian model.
        ///
        /// EFD inverse for the ordinary Hessian:
        ///
        ///   $$x' = (36*(delta^3 - 1) - v)/(6*(u + 9*delta^2)) - delta/2$$
        ///   $$y' = (v + 36*(delta^3 - 1))/(6*(u + 9*delta^2)) - delta/2$$
        ///
        /// Then recover the generalized-Hessian coordinates by x = zeta*x', y = zeta*y'.
        pub fn map_point_from_weierstrass_with_zeta(
            &self,
            p: &AffinePoint<F>,
            zeta: F,
        ) -> Option<HessianPoint<F>> {
            if p.infinity {
                return Some(HessianPoint::identity());
            }

            let delta = self.efd_parameter_with_zeta(zeta)?;
            let delta_square = <F as FieldOps>::square(&delta);
            let delta_cubic = &delta * &delta_square;
            let dm1 = &delta_cubic - &F::one();

            let two_inv =  F::from_u64(2).invert().into_option()?;
            let denom = &F::from_u64(6) * &(&p.x + &(&F::from_u64(9) * &delta_square));
            let denom_inv = match denom.invert().into_option() {
                Some(v) => v,
                None => return Some(HessianPoint::identity()),
            };

            let common = &F::from_u64(36) * &dm1;
            let xh = &(&(&common - &p.y) * &denom_inv) - &(&delta * &two_inv);
            let yh = &(&(&p.y + &common) * &denom_inv) - &(&delta * &two_inv);

            Some(HessianPoint::from_affine(&zeta * &xh, &zeta * &yh))
        }
    }
}
