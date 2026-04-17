//! Tests for Legendre point operations (group law and scalar multiplication).
//!
//! We work over F_19 with the Legendre curve
//!
//!     y^2 = x(x-1)(x-3).
//!
//! This curve has 16 rational points in total, and the point P = (2,6)
//! has order 8. The visible 2-torsion points are (0,0), (1,0), and (3,0).

use crypto_bigint::{const_prime_monty_params, Uint};

use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

use ec::curve_legendre::LegendreCurve;
use ec::curve_ops::Curve;
use ec::point_legendre::LegendrePoint;

// ===========================================================================
// Field definitions
// ===========================================================================

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(n: u64) -> F19 {
    F19::from_u64(n)
}

fn curve() -> LegendreCurve<F19> {
    LegendreCurve::new(fp(3))
}

// ===========================================================================
// Helpers
// ===========================================================================

/// Brute-force all affine points on y^2 = x(x-1)(x-3) over F_19.
fn all_legendre_19(c: &LegendreCurve<F19>) -> Vec<LegendrePoint<F19>> {
    let mut pts = Vec::new();
    for xv in 0..19u64 {
        for yv in 0..19u64 {
            let p = LegendrePoint::new(fp(xv), fp(yv));
            if c.is_on_curve(&p) {
                pts.push(p);
            }
        }
    }
    pts
}

/// A point of order 8 on the test curve.
fn p() -> LegendrePoint<F19> {
    LegendrePoint::new(fp(2), fp(6))
}

// ===========================================================================
// Point tests over F_19
// ===========================================================================

#[test]
fn f19_identity_is_identity() {
    let id = LegendrePoint::<F19>::identity();

    assert!(id.is_identity());
    assert!(id.infinity);
    assert_eq!(id.x, F19::zero());
    assert_eq!(id.y, F19::one());
}
#[test]
fn f19_add_identity() {
    let c = curve();
    let id = LegendrePoint::<F19>::identity();
    let p = p();

    assert_eq!(p.add(&id, &c), p);
    assert_eq!(id.add(&p, &c), p);
}

#[test]
fn f19_negation() {
    let c = curve();
    let p = p();
    let neg = p.negate(&c);

    assert_eq!(neg, LegendrePoint::new(fp(2), fp(13)));
    assert!(p.add(&neg, &c).is_identity());
}

#[test]
fn f19_negate_then_add_for_all_points() {
    let c = curve();
    for p in all_legendre_19(&c) {
        let neg = p.negate(&c);
        assert!(
            p.add(&neg, &c).is_identity(),
            "P + (-P) != O for P = {}",
            p
        );
    }
}

#[test]
fn f19_visible_two_torsion_doubles_to_identity() {
    let c = curve();

    let t0 = LegendrePoint::new(fp(0), fp(0));
    let t1 = LegendrePoint::new(fp(1), fp(0));
    let tl = LegendrePoint::new(fp(3), fp(0));

    assert_eq!(t0.double(&c), LegendrePoint::identity());
    assert_eq!(t1.double(&c), LegendrePoint::identity());
    assert_eq!(tl.double(&c), LegendrePoint::identity());
}

#[test]
fn f19_double_known_vector() {
    let c = curve();
    let p = p();

    // 2*(2,6) = (7,15)
    let p2 = p.double(&c);
    assert_eq!(p2, LegendrePoint::new(fp(7), fp(15)));
    assert!(c.is_on_curve(&p2));
}

#[test]
fn f19_add_known_vector() {
    let c = curve();
    let p = p();
    let q = LegendrePoint::new(fp(7), fp(15));

    // (2,6) + (7,15) = (18,7)
    let sum = p.add(&q, &c);
    assert_eq!(sum, LegendrePoint::new(fp(18), fp(7)));
    assert!(c.is_on_curve(&sum));
}

#[test]
fn f19_add_self_matches_double() {
    let c = curve();
    let p = p();

    assert_eq!(p.add(&p, &c), p.double(&c));
}

#[test]
fn f19_double_on_curve_for_all_points() {
    let c = curve();
    for p in all_legendre_19(&c) {
        let dbl = p.double(&c);
        assert!(c.is_on_curve(&dbl));
    }
}

#[test]
fn f19_add_on_curve_small_sample() {
    let c = curve();
    let pts = all_legendre_19(&c);

    for i in 0..pts.len().min(8) {
        for j in 0..pts.len().min(8) {
            let s = pts[i].add(&pts[j], &c);
            assert!(c.is_on_curve(&s));
        }
    }
}

#[test]
fn f19_commutativity_small_sample() {
    let c = curve();
    let pts = all_legendre_19(&c);

    for i in 0..pts.len().min(8) {
        for j in (i + 1)..pts.len().min(8) {
            assert_eq!(pts[i].add(&pts[j], &c), pts[j].add(&pts[i], &c));
        }
    }
}

#[test]
fn f19_associativity_small_sample() {
    let c = curve();
    let pts = all_legendre_19(&c);

    assert!(pts.len() >= 3);
    let p = pts[0];
    let q = pts[1];
    let r = pts[2];

    assert_eq!(
        p.add(&q, &c).add(&r, &c),
        p.add(&q.add(&r, &c), &c),
    );
}

#[test]
fn f19_scalar_mul_order() {
    let c = curve();
    let p = p();

    // P = (2,6) has order 8 on this curve.
    assert!(p.scalar_mul(&[8], &c).is_identity());
}

#[test]
fn f19_scalar_mul_consistency() {
    let c = curve();
    let p = p();

    let seven_p = p.scalar_mul(&[7], &c);
    let mut acc = LegendrePoint::<F19>::identity();
    for _ in 0..7 {
        acc = acc.add(&p, &c);
    }

    assert_eq!(seven_p, acc);
}

#[test]
fn f19_scalar_mul_double() {
    let c = curve();
    let p = p();

    assert_eq!(p.scalar_mul(&[2], &c), p.double(&c));
}

#[test]
fn f19_chain_on_curve() {
    let c = curve();
    let p = p();

    let mut acc = p;
    for _ in 1..16 {
        acc = acc.add(&p, &c);
        assert!(c.is_on_curve(&acc));
    }
}