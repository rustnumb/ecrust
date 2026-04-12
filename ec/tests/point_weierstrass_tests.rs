//! Integration tests for the `ec` crate.
//!
//! Tests cover:
//!   1. Short Weierstrass over a prime field  (F₁₉)
//!   2. Short Weierstrass over a quadratic extension  (F₁₉²)
//!   3. General Weierstrass over a binary extension field  (F₂⁴)

use crypto_bigint::{const_prime_monty_params, Uint};

use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;
use fp::fp_ext::{FpExt, IrreduciblePoly, TonelliShanksConstants};

use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_weierstrass::AffinePoint;

// ===========================================================================
// 1.  Short Weierstrass over F₁₉
// ===========================================================================

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;
fn fp(n: u64) -> F19 {
    F19::from_u64(n)
}

/// Brute-force: all affine points on  y² = x³ + ax + b  over F_p.
fn all_points_short_19(a: u64, b: u64) -> Vec<(u64, u64)> {
    let p = 19u64;
    let mut pts = Vec::new();
    for x in 0..p {
        let rhs = (x * x % p * x % p + a * x % p + b) % p;
        for y in 0..p {
            if (y * y) % p == rhs {
                pts.push((x, y));
            }
        }
    }
    pts
}

#[test]
fn fp19_identity_is_at_infinity() {
    assert!(AffinePoint::<F19>::identity().is_identity());
}

#[test]
fn fp19_negate_then_add_is_identity() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    for (xv, yv) in all_points_short_19(2, 3) {
        let p = AffinePoint::new(fp(xv), fp(yv));
        let neg = p.negate(&c);
        assert!(
            p.add(&neg, &c).is_identity(),
            "P + (−P) should be O for ({}, {})",
            xv,
            yv,
        );
    }
}

#[test]
fn fp19_commutativity() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    let pts = all_points_short_19(2, 3);
    for i in 0..pts.len().min(5) {
        for j in (i + 1)..pts.len().min(5) {
            let p = AffinePoint::new(fp(pts[i].0), fp(pts[i].1));
            let q = AffinePoint::new(fp(pts[j].0), fp(pts[j].1));
            assert_eq!(
                p.add(&q, &c),
                q.add(&p, &c),
                "P+Q ≠ Q+P for P=({},{}), Q=({},{})",
                pts[i].0,
                pts[i].1,
                pts[j].0,
                pts[j].1,
            );
        }
    }
}

#[test]
fn fp19_associativity() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    let pts = all_points_short_19(2, 3);
    assert!(pts.len() >= 3);
    let p = AffinePoint::new(fp(pts[0].0), fp(pts[0].1));
    let q = AffinePoint::new(fp(pts[1].0), fp(pts[1].1));
    let r = AffinePoint::new(fp(pts[2].0), fp(pts[2].1));
    assert_eq!(p.add(&q, &c).add(&r, &c), p.add(&q.add(&r, &c), &c),);
}

#[test]
fn fp19_scalar_mul_order() {
    let pts = all_points_short_19(2, 3);
    let order = pts.len() as u64 + 1;

    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    let y = fp(6).sqrt().into_option().expect("6 must be a QR mod 19");
    let p = AffinePoint::new(fp(1), y);
    assert!(
        p.scalar_mul(&[order], &c).is_identity(),
        "[#E]P must be O  (order = {})",
        order,
    );
}

#[test]
fn fp19_scalar_mul_consistency() {
    // [5]P should equal P+P+P+P+P
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    let y = fp(6).sqrt().into_option().expect("6 must be a QR mod 19");
    let p = AffinePoint::new(fp(1), y);

    let five_p = p.scalar_mul(&[5], &c);
    let mut acc = AffinePoint::identity();
    for _ in 0..5 {
        acc = acc.add(&p, &c);
    }
    assert_eq!(five_p, acc);
}

// ===========================================================================
// 2.  Short Weierstrass over F₁₉²  (extension field)
// ===========================================================================

struct QuadPoly;
struct TSQuad;

impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
    fn modulus() -> [F19; 2] {
        [F19::one(), F19::zero()] // x² + 1
    }
}

impl TonelliShanksConstants<Fp19Mod, 1, 2, 1> for TSQuad {
    // Still only need 1 limb for 19^2
    const ORDER: Uint<1> = Uint::<1>::from_u64(360);
    const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(180);
    const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(22);
    const TWOSM1: Uint<1> = Uint::<1>::from_u64(4);
    fn root_of_unity() -> [FpElement<Fp19Mod, 1>; 2] {
        [F19::from_u64(3), F19::from_u64(3)]
    }
    const S: u64 = 3;
    const T: Uint<1> = Uint::<1>::from_u64(45);
}

type F19_2 = FpExt<Fp19Mod, 1, 2, 1, QuadPoly, TSQuad>;

fn el(a: u64, b: u64) -> F19_2 {
    F19_2::new([fp(a), fp(b)])
}

/// Try to find a point on y² = x³ + ax + b over F₁₉² by scanning a few
/// elements and testing whether the RHS is a square (via Euler criterion).
fn find_point_fp19_2(curve: &WeierstrassCurve<F19_2>) -> Option<AffinePoint<F19_2>> {
    // |F₁₉²*| = 19² − 1 = 360,  so  (p^m−1)/2 = 180
    let half_order: &[u64] = &[180];

    for a0 in 0..19u64 {
        for a1 in 0..19u64 {
            let x = el(a0, a1);
            let x2 = <F19_2 as FieldOps>::square(&x);
            let x3 = <F19_2 as FieldOps>::mul(&x2, &x);
            let rhs = <F19_2 as FieldOps>::add(
                &<F19_2 as FieldOps>::add(&x3, &<F19_2 as FieldOps>::mul(&curve.a4, &x)),
                &curve.a6,
            );
            if bool::from(rhs.is_zero()) {
                return Some(AffinePoint::new(x, F19_2::zero()));
            }
            // Euler criterion: rhs is a QR iff rhs^((q-1)/2) == 1
            let euler = rhs.pow(half_order);
            if bool::from(euler.is_one()) {
                // We know a square root exists.  Find it by brute force
                // over the small field.
                for b0 in 0..19u64 {
                    for b1 in 0..19u64 {
                        let y = el(b0, b1);
                        if <F19_2 as FieldOps>::square(&y) == rhs {
                            return Some(AffinePoint::new(x, y));
                        }
                    }
                }
            }
        }
    }
    None
}

#[test]
fn fp19_2_point_on_curve() {
    // y² = x³ + (1+0x)·x + (2+0x)  over F₁₉²
    let c = WeierstrassCurve::new_short(el(1, 0), el(2, 0));
    let p = find_point_fp19_2(&c).expect("should find at least one point");
    assert!(c.contains(&p.x, &p.y));
}

#[test]
fn fp19_2_double_on_curve() {
    let c = WeierstrassCurve::new_short(el(1, 0), el(2, 0));
    let p = find_point_fp19_2(&c).expect("need a point");
    let dbl = p.double(&c);
    if !dbl.is_identity() {
        assert!(c.contains(&dbl.x, &dbl.y));
    }
}

#[test]
fn fp19_2_add_inverse() {
    let c = WeierstrassCurve::new_short(el(1, 0), el(2, 0));
    let p = find_point_fp19_2(&c).expect("need a point");
    let neg = p.negate(&c);
    assert!(p.add(&neg, &c).is_identity());
}

// ===========================================================================
// 3.  General Weierstrass over F₂⁴  (binary extension)
// ===========================================================================

use fp::f2_ext::{BinaryIrreducible, F2Ext};

/// Irreducible polynomial  x⁴ + x + 1  over F₂.
/// Binary representation: 10011 = 0x13.
struct Poly2_4;
impl BinaryIrreducible<1> for Poly2_4 {
    fn modulus() -> Uint<1> {
        Uint::from_u64(0x13)
    } // x⁴ + x + 1
    fn degree() -> usize {
        4
    }
}

type GF16 = F2Ext<1, Poly2_4>;

fn gf(n: u64) -> GF16 {
    GF16::from_u64(n)
}

/// In characteristic 2 the short Weierstrass form degenerates, so we use
/// the general Weierstrass form:
///
/// ```text
/// y² + xy = x³ + a₂x² + a₆        (a₁=1, a₃=0, a₄=0)
/// ```
///
/// This is a standard "non-supersingular" form for binary curves.
fn binary_curve() -> WeierstrassCurve<GF16> {
    WeierstrassCurve::new(
        gf(1),      // a₁ = 1
        gf(0b1000), // a₂ = x^3   (some nonzero element)
        gf(0),      // a₃ = 0
        gf(0),      // a₄ = 0
        gf(1),      // a₆ = 1
    )
}

/// Brute-force search for all affine points on the binary curve over GF(16).
fn all_points_gf16(c: &WeierstrassCurve<GF16>) -> Vec<AffinePoint<GF16>> {
    let mut pts = Vec::new();
    for xv in 0..16u64 {
        for yv in 0..16u64 {
            let x = gf(xv);
            let y = gf(yv);
            if c.contains(&x, &y) {
                pts.push(AffinePoint::new(x, y));
            }
        }
    }
    pts
}

#[test]
fn gf16_find_points() {
    let c = binary_curve();
    let pts = all_points_gf16(&c);
    assert!(
        !pts.is_empty(),
        "expected at least one affine point on binary curve",
    );
}

#[test]
fn gf16_double_on_curve() {
    let c = binary_curve();
    for p in all_points_gf16(&c) {
        let dbl = p.double(&c);
        if !dbl.is_identity() {
            assert!(
                c.contains(&dbl.x, &dbl.y),
                "[2]P not on curve for P=({:?}, {:?})",
                p.x.as_uint().to_words()[0],
                p.y.as_uint().to_words()[0],
            );
        }
    }
}

#[test]
fn gf16_add_inverse() {
    let c = binary_curve();
    for p in all_points_gf16(&c) {
        let neg = p.negate(&c);
        assert!(
            p.add(&neg, &c).is_identity(),
            "P + (−P) ≠ O for P=({:?}, {:?})",
            p.x.as_uint().to_words()[0],
            p.y.as_uint().to_words()[0],
        );
    }
}

#[test]
fn gf16_add_on_curve() {
    let c = binary_curve();
    let pts = all_points_gf16(&c);
    for i in 0..pts.len().min(8) {
        for j in 0..pts.len().min(8) {
            let r = pts[i].add(&pts[j], &c);
            if !r.is_identity() {
                assert!(c.contains(&r.x, &r.y), "P+Q not on curve",);
            }
        }
    }
}

#[test]
fn gf16_commutativity() {
    let c = binary_curve();
    let pts = all_points_gf16(&c);
    for i in 0..pts.len().min(6) {
        for j in (i + 1)..pts.len().min(6) {
            assert_eq!(pts[i].add(&pts[j], &c), pts[j].add(&pts[i], &c),);
        }
    }
}

#[test]
fn gf16_associativity() {
    let c = binary_curve();
    let pts = all_points_gf16(&c);
    if pts.len() >= 3 {
        let p = &pts[0];
        let q = &pts[1];
        let r = &pts[2];
        assert_eq!(p.add(q, &c).add(r, &c), p.add(&q.add(r, &c), &c),);
    }
}

#[test]
fn gf16_scalar_mul_order() {
    let c = binary_curve();
    let pts = all_points_gf16(&c);
    let group_order = pts.len() as u64 + 1; // affine points + O

    for p in &pts {
        assert!(
            p.scalar_mul(&[group_order], &c).is_identity(),
            "[{}]P should be O",
            group_order,
        );
    }
}