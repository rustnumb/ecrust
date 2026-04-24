//! Tests for Edwards curve construction, membership, and invariants.
//!
//! Covers:
//!   1. F₇₆₁   (p = 761,   d = 3, non-square → complete)
//!   2. F₆₅₅₃₇ (p = 65537, d = 3, non-square → complete)
//!   3. GF(2⁴)  (d₁ = 1, d₂ = α³ = 8, Tr(d₂) = 1 → complete)
//!
//! Sage verification sketch:
//!
//! ```sage
//! # F_761
//! p=761; F=GF(p); assert not F(3).is_square()
//! # F_65537
//! p=65537; F=GF(p); assert not F(3).is_square()
//! # GF(2^4)
//! R.<x>=GF(2)[]; F.<a>=GF(2^4,modulus=x^4+x+1); assert (a^3).trace()==1
//! ```

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
// 1.  F₇₆₁, d = 3
// ===========================================================================

#[test]
fn f761_construction() {
    let c = EdwardsCurve::new(fp761(3));
    assert_eq!(c.d(), fp761(3));
}

#[test]
#[should_panic]
fn f761_d_zero_panics() {
    let _ = EdwardsCurve::new(fp761(0));
}

#[test]
#[should_panic]
fn f761_d_one_panics() {
    let _ = EdwardsCurve::new(fp761(1));
}

#[test]
fn f761_identity_on_curve() {
    let c = EdwardsCurve::new(fp761(3));
    let id = EdwardsPoint::<F761>::identity();
    assert!(c.is_on_curve(&id));
}

#[test]
fn f761_known_points_on_curve() {
    let c = EdwardsCurve::new(fp761(3));

    // (1,0): 1+0 = 1+0 ✓
    assert!(c.contains(&fp761(1), &fp761(0)));
    // (760,0) = (-1,0): 1+0 = 1+0 ✓
    assert!(c.contains(&fp761(760), &fp761(0)));
    // (0,760) = (0,-1): 0+1 = 1+0 ✓
    assert!(c.contains(&fp761(0), &fp761(760)));
}

#[test]
fn f761_rejects_non_point() {
    let c = EdwardsCurve::new(fp761(3));
    // (1,1): 1+1=2, 1+3=4, 2≠4
    assert!(!c.contains(&fp761(1), &fp761(1)));
}

#[test]
fn f761_j_invariant_nonzero() {
    let c = EdwardsCurve::new(fp761(3));
    assert!(!bool::from(c.j_invariant().is_zero()));
}

#[test]
fn f761_a_invariants() {
    let c = EdwardsCurve::new(fp761(3));
    let inv = c.a_invariants();
    assert_eq!(inv.len(), 1);
    assert_eq!(inv[0], fp761(3));
}

// ===========================================================================
// 2.  F₆₅₅₃₇, d = 3
// ===========================================================================

#[test]
fn f65537_construction() {
    let c = EdwardsCurve::new(fp65537(3));
    assert_eq!(c.d(), fp65537(3));
}

#[test]
fn f65537_identity_on_curve() {
    let c = EdwardsCurve::new(fp65537(3));
    assert!(c.is_on_curve(&EdwardsPoint::<F65537>::identity()));
}

#[test]
fn f65537_known_points_on_curve() {
    let c = EdwardsCurve::new(fp65537(3));
    assert!(c.contains(&fp65537(1), &fp65537(0)));
    assert!(c.contains(&fp65537(65536), &fp65537(0)));
    assert!(c.contains(&fp65537(0), &fp65537(65536)));
}

#[test]
fn f65537_rejects_non_point() {
    let c = EdwardsCurve::new(fp65537(3));
    assert!(!c.contains(&fp65537(2), &fp65537(2)));
}

#[test]
fn f65537_j_invariant_nonzero() {
    let c = EdwardsCurve::new(fp65537(3));
    assert!(!bool::from(c.j_invariant().is_zero()));
}

// ===========================================================================
// 3.  Binary Edwards over GF(2⁴),  d₁ = 1, d₂ = 8
// ===========================================================================

fn binary_edwards_curve() -> EdwardsCurve<GF16> {
    EdwardsCurve::new_binary(gf(1), gf(8))
}

#[test]
fn gf16_construction() {
    let c = binary_edwards_curve();
    assert_eq!(c.d1, gf(1));
    assert_eq!(c.d2, gf(8));
}

#[test]
#[should_panic]
fn gf16_d1_zero_panics() {
    let _ = EdwardsCurve::new_binary(gf(0), gf(8));
}

#[test]
fn gf16_identity_on_curve() {
    let c = binary_edwards_curve();
    let id = EdwardsPoint::<GF16>::identity();
    assert!(c.is_on_curve(&id));
    assert_eq!(id, EdwardsPoint::new(gf(0), gf(0)));
}

#[test]
fn gf16_find_points() {
    let c = binary_edwards_curve();
    let mut count = 0usize;
    for xv in 0..16u64 {
        for yv in 0..16u64 {
            if c.contains(&gf(xv), &gf(yv)) {
                count += 1;
            }
        }
    }
    assert!(count > 0, "expected at least one point");
}

#[test]
fn gf16_order2_point_on_curve() {
    // (1,1) should lie on the curve and have order 2
    let c = binary_edwards_curve();
    assert!(c.contains(&gf(1), &gf(1)));
}

#[test]
fn gf16_j_invariant_nonzero() {
    let c = binary_edwards_curve();
    assert!(!bool::from(c.j_invariant().is_zero()));
}

#[test]
fn gf16_a_invariants() {
    let c = binary_edwards_curve();
    let inv = c.a_invariants();
    assert_eq!(inv.len(), 2);
    assert_eq!(inv[0], gf(1));
    assert_eq!(inv[1], gf(8));
}
