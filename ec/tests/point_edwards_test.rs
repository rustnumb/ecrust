//! Tests for Edwards point operations (group law, scalar mul, w-coords).
//!
//! Covers:
//!   1. F₇₆₁   (p = 761,   d = 3, non-square → complete)
//!   2. F₆₅₅₃₇ (p = 65537, d = 3, non-square → complete)
//!   3. GF(2⁴)  (d₁ = 1, d₂ = α³ = 8, Tr(d₂) = 1 → complete)

use crypto_bigint::{Uint, const_prime_monty_params};

use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

use ec::curve_edwards::EdwardsCurve;
use ec::curve_ops::Curve;
use ec::point_edwards::EdwardsPoint;

// ===========================================================================
// Field definitions
// ===========================================================================

const_prime_monty_params!(Fp761Mod, Uint<1>, "00000000000002F9", 2);
type F761 = FpElement<Fp761Mod, 1>;
fn fp761(n: u64) -> F761 {
    F761::from_u64(n)
}

const_prime_monty_params!(Fp65537Mod, Uint<1>, "0000000000010001", 3);
type F65537 = FpElement<Fp65537Mod, 1>;
fn fp65537(n: u64) -> F65537 {
    F65537::from_u64(n)
}

use fp::f2_ext::{BinaryIrreducible, F2Ext};

struct Poly2_4;
impl BinaryIrreducible<1> for Poly2_4 {
    fn modulus() -> Uint<1> {
        Uint::from_u64(0x13)
    }
    fn degree() -> usize {
        4
    }
}
type GF16 = F2Ext<1, Poly2_4>;
fn gf(n: u64) -> GF16 {
    GF16::from_u64(n)
}

// ===========================================================================
// Helpers
// ===========================================================================

/// Brute-force all affine points on x²+y²=1+d·x²·y² over F₇₆₁.
fn all_edwards_761(d: u64) -> Vec<(u64, u64)> {
    let p = 761u64;
    let mut pts = Vec::new();
    for x in 0..p {
        for y in 0..p {
            let x2 = x * x % p;
            let y2 = y * y % p;
            let lhs = (x2 + y2) % p;
            let rhs = (1 + d % p * x2 % p * y2 % p) % p;
            if lhs == rhs {
                pts.push((x, y));
            }
        }
    }
    pts
}

/// Find a second point on Edwards over F₆₅₅₃₇ by scanning x.
fn find_point_65537(curve: &EdwardsCurve<F65537>) -> Option<EdwardsPoint<F65537>> {
    let d = curve.d();
    for xv in 2u64..1000 {
        let x = fp65537(xv);
        let x2 = <F65537 as FieldOps>::square(&x);
        let num = &F65537::one() - &x2;
        let den = &F65537::one() - &(&d * &x2);
        if let Some(den_inv) = den.invert().into_option() {
            let y2 = &num * &den_inv;
            if let Some(y) = y2.sqrt().into_option() {
                let p = EdwardsPoint::new(x, y);
                if curve.is_on_curve(&p) {
                    return Some(p);
                }
            }
        }
    }
    None
}

/// Brute-force all affine points on binary Edwards over GF(16).
fn all_bin_edwards(c: &EdwardsCurve<GF16>) -> Vec<EdwardsPoint<GF16>> {
    let mut pts = Vec::new();
    for xv in 0..16u64 {
        for yv in 0..16u64 {
            if c.contains(&gf(xv), &gf(yv)) {
                pts.push(EdwardsPoint::new(gf(xv), gf(yv)));
            }
        }
    }
    pts
}

fn binary_edwards_curve() -> EdwardsCurve<GF16> {
    EdwardsCurve::new_binary(gf(1), gf(8))
}

// ===========================================================================
// 1.  Point tests over F₇₆₁
// ===========================================================================

#[test]
fn f761_identity_is_identity() {
    assert!(EdwardsPoint::<F761>::identity().is_identity());
}

#[test]
fn f761_add_identity() {
    let c = EdwardsCurve::new(fp761(3));
    let id = EdwardsPoint::<F761>::identity();
    let p = EdwardsPoint::new(fp761(1), fp761(0));
    assert_eq!(p.add(&id, &c), p);
    assert_eq!(id.add(&p, &c), p);
}

#[test]
fn f761_negation() {
    let c = EdwardsCurve::new(fp761(3));
    let p = EdwardsPoint::new(fp761(1), fp761(0));
    let neg = p.negate(&c);
    // -(1,0) = (-1,0) = (760,0)
    assert_eq!(neg, EdwardsPoint::new(fp761(760), fp761(0)));
}

#[test]
fn f761_negate_then_add() {
    let c = EdwardsCurve::new(fp761(3));
    for (xv, yv) in all_edwards_761(3) {
        let p = EdwardsPoint::new(fp761(xv), fp761(yv));
        let neg = p.negate(&c);
        assert!(
            p.add(&neg, &c).is_identity(),
            "P + (-P) ≠ O for ({},{})",
            xv,
            yv,
        );
    }
}

#[test]
fn f761_doubling_order4() {
    // (1,0) has order 4: 2(1,0)=(0,760), 4(1,0)=O
    let c = EdwardsCurve::new(fp761(3));
    let p = EdwardsPoint::new(fp761(1), fp761(0));

    let p2 = p.double(&c);
    assert_eq!(p2, EdwardsPoint::new(fp761(0), fp761(760)));

    let p4 = p2.double(&c);
    assert!(p4.is_identity());
}

#[test]
fn f761_double_on_curve() {
    let c = EdwardsCurve::new(fp761(3));
    for (xv, yv) in all_edwards_761(3) {
        let p = EdwardsPoint::new(fp761(xv), fp761(yv));
        let dbl = p.double(&c);
        assert!(c.is_on_curve(&dbl));
    }
}

#[test]
fn f761_commutativity() {
    let c = EdwardsCurve::new(fp761(3));
    let pts = all_edwards_761(3);
    for i in 0..pts.len().min(6) {
        for j in (i + 1)..pts.len().min(6) {
            let p = EdwardsPoint::new(fp761(pts[i].0), fp761(pts[i].1));
            let q = EdwardsPoint::new(fp761(pts[j].0), fp761(pts[j].1));
            assert_eq!(p.add(&q, &c), q.add(&p, &c));
        }
    }
}

#[test]
fn f761_associativity() {
    let c = EdwardsCurve::new(fp761(3));
    let pts = all_edwards_761(3);
    assert!(pts.len() >= 3);
    let p = EdwardsPoint::new(fp761(pts[0].0), fp761(pts[0].1));
    let q = EdwardsPoint::new(fp761(pts[1].0), fp761(pts[1].1));
    let r = EdwardsPoint::new(fp761(pts[2].0), fp761(pts[2].1));
    assert_eq!(p.add(&q, &c).add(&r, &c), p.add(&q.add(&r, &c), &c),);
}

#[test]
fn f761_scalar_mul_order() {
    let c = EdwardsCurve::new(fp761(3));
    let pts = all_edwards_761(3);
    let order = pts.len() as u64;
    let p = EdwardsPoint::new(fp761(pts[0].0), fp761(pts[0].1));
    assert!(
        p.scalar_mul(&[order], &c).is_identity(),
        "[{}]P should be O",
        order,
    );
}

#[test]
fn f761_scalar_mul_consistency() {
    let c = EdwardsCurve::new(fp761(3));
    let p = EdwardsPoint::new(fp761(1), fp761(0));
    let five_p = p.scalar_mul(&[5], &c);
    let mut acc = EdwardsPoint::<F761>::identity();
    for _ in 0..5 {
        acc = acc.add(&p, &c);
    }
    assert_eq!(five_p, acc);
}

#[test]
fn f761_scalar_mul_double() {
    let c = EdwardsCurve::new(fp761(3));
    let p = EdwardsPoint::new(fp761(1), fp761(0));
    assert_eq!(p.scalar_mul(&[2], &c), p.double(&c));
}

// ===========================================================================
// 2.  Point tests over F₆₅₅₃₇
// ===========================================================================

#[test]
fn f65537_identity_is_identity() {
    assert!(EdwardsPoint::<F65537>::identity().is_identity());
    assert_eq!(
        EdwardsPoint::<F65537>::identity(),
        EdwardsPoint::new(fp65537(0), fp65537(1)),
    );
}

#[test]
fn f65537_add_identity() {
    let c = EdwardsCurve::new(fp65537(3));
    let id = EdwardsPoint::<F65537>::identity();
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    assert_eq!(p.add(&id, &c), p);
    assert_eq!(id.add(&p, &c), p);
}

#[test]
fn f65537_negation() {
    let c = EdwardsCurve::new(fp65537(3));
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    let neg = p.negate(&c);
    assert_eq!(neg, EdwardsPoint::new(fp65537(65536), fp65537(0)));
    assert!(p.add(&neg, &c).is_identity());
}

#[test]
fn f65537_order4() {
    let c = EdwardsCurve::new(fp65537(3));
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    let p2 = p.double(&c);
    assert_eq!(p2, EdwardsPoint::new(fp65537(0), fp65537(65536)));
    assert!(p2.double(&c).is_identity());
}

#[test]
fn f65537_commutativity() {
    let c = EdwardsCurve::new(fp65537(3));
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    let q = match find_point_65537(&c) {
        Some(q) => q,
        None => p.double(&c),
    };
    assert_eq!(p.add(&q, &c), q.add(&p, &c));
}

#[test]
fn f65537_associativity() {
    let c = EdwardsCurve::new(fp65537(3));
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    let q = p.double(&c);
    let r = q.double(&c);
    assert_eq!(p.add(&q, &c).add(&r, &c), p.add(&q.add(&r, &c), &c),);
}

#[test]
fn f65537_scalar_mul_consistency() {
    let c = EdwardsCurve::new(fp65537(3));
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    let seven_p = p.scalar_mul(&[7], &c);
    let mut acc = EdwardsPoint::<F65537>::identity();
    for _ in 0..7 {
        acc = acc.add(&p, &c);
    }
    assert_eq!(seven_p, acc);
}

#[test]
fn f65537_chain_on_curve() {
    let c = EdwardsCurve::new(fp65537(3));
    let p = EdwardsPoint::new(fp65537(1), fp65537(0));
    let mut acc = p;
    for _ in 1..20 {
        acc = acc.add(&p, &c);
        assert!(c.is_on_curve(&acc));
    }
}

// ===========================================================================
// 3.  Point tests over GF(2⁴), binary Edwards
// ===========================================================================

#[test]
fn gf16_identity_is_identity() {
    let id = EdwardsPoint::<GF16>::identity();
    assert!(id.is_identity());
    assert_eq!(id, EdwardsPoint::new(gf(0), gf(0)));
}

#[test]
fn gf16_negation_swaps() {
    let c = binary_edwards_curve();
    for p in all_bin_edwards(&c) {
        let neg = p.negate(&c);
        assert_eq!(neg.x, p.y);
        assert_eq!(neg.y, p.x);
    }
}

#[test]
fn gf16_negate_then_add() {
    let c = binary_edwards_curve();
    for p in all_bin_edwards(&c) {
        assert!(
            p.add(&p.negate(&c), &c).is_identity(),
            "P+(-P) ≠ O for ({:?},{:?})",
            p.x.as_uint().to_words()[0],
            p.y.as_uint().to_words()[0],
        );
    }
}

#[test]
fn gf16_add_identity() {
    let c = binary_edwards_curve();
    let id = EdwardsPoint::<GF16>::identity();
    for p in all_bin_edwards(&c) {
        assert_eq!(p.add(&id, &c), p);
        assert_eq!(id.add(&p, &c), p);
    }
}

#[test]
fn gf16_order2_point() {
    let c = binary_edwards_curve();
    let p = EdwardsPoint::new(gf(1), gf(1));
    assert!(c.is_on_curve(&p));
    assert_eq!(p.negate(&c), p, "(1,1) is its own negation");
    assert!(p.double(&c).is_identity(), "2*(1,1) should be O");
}

#[test]
fn gf16_double_on_curve() {
    let c = binary_edwards_curve();
    for p in all_bin_edwards(&c) {
        assert!(c.is_on_curve(&p.double(&c)));
    }
}

#[test]
fn gf16_add_on_curve() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    for i in 0..pts.len().min(8) {
        for j in 0..pts.len().min(8) {
            assert!(c.is_on_curve(&pts[i].add(&pts[j], &c)));
        }
    }
}

#[test]
fn gf16_commutativity() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    for i in 0..pts.len().min(8) {
        for j in (i + 1)..pts.len().min(8) {
            assert_eq!(pts[i].add(&pts[j], &c), pts[j].add(&pts[i], &c));
        }
    }
}

#[test]
fn gf16_associativity() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    if pts.len() >= 3 {
        assert_eq!(
            pts[0].add(&pts[1], &c).add(&pts[2], &c),
            pts[0].add(&pts[1].add(&pts[2], &c), &c),
        );
    }
}

#[test]
fn gf16_scalar_mul_order() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    let order = pts.len() as u64;
    for p in &pts {
        assert!(
            p.scalar_mul(&[order], &c).is_identity(),
            "[{}]P should be O",
            order,
        );
    }
}

#[test]
fn gf16_scalar_mul_consistency() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    if pts.len() >= 2 {
        let p = &pts[1];
        let five_p = p.scalar_mul(&[5], &c);
        let mut acc = EdwardsPoint::<GF16>::identity();
        for _ in 0..5 {
            acc = acc.add(p, &c);
        }
        assert_eq!(five_p, acc);
    }
}

#[test]
fn gf16_scalar_mul_double() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    for p in &pts {
        assert_eq!(p.scalar_mul(&[2], &c), p.double(&c));
    }
}

#[test]
fn gf16_w_double() {
    let c = binary_edwards_curve();
    for p in all_bin_edwards(&c) {
        let w = &p.x + &p.y;
        let dbl = p.double(&c);
        let w_expected = &dbl.x + &dbl.y;
        let w_got = EdwardsPoint::<GF16>::w_double(&w, &c);
        assert_eq!(
            w_got,
            w_expected,
            "w-double mismatch for P=({:?},{:?})",
            p.x.as_uint().to_words()[0],
            p.y.as_uint().to_words()[0],
        );
    }
}

#[test]
fn gf16_w_diff_add() {
    let c = binary_edwards_curve();
    let pts = all_bin_edwards(&c);
    if pts.len() >= 3 {
        // Pick P, Q such that Q-P is also known
        let p = &pts[1];
        let q = &pts[2];
        let diff = q.add(&p.negate(&c), &c); // Q - P
        let sum = p.add(q, &c); // P + Q

        let w1 = &diff.x + &diff.y;
        let w2 = &p.x + &p.y;
        let w3 = &q.x + &q.y;
        let w5_expected = &sum.x + &sum.y;

        let w5_got = EdwardsPoint::<GF16>::w_diff_add(&w1, &w2, &w3, &c);
        assert_eq!(w5_got, w5_expected);
    }
}
