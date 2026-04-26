//! Integration tests for twisted Hessian curves.

use crypto_bigint::{Uint, const_prime_monty_params};

use ec::curve_ops::Curve;
use ec::curve_twisted_hessian::TwistedHessianCurve;
use ec::point_ops::{PointAdd, PointOps};
use ec::point_twisted_hessian::TwistedHessianPoint;
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp103Mod, Uint<1>, "0000000000000067", 2);
type F103 = FpElement<Fp103Mod, 1>;

fn f(n: i64) -> F103 {
    F103::from_u64(n.rem_euclid(103) as u64)
}

fn test_curve() -> TwistedHessianCurve<F103> {
    // Twisted normal form a = 8, d = 1.
    // Since 1^3 != 8 in F_103, this curve is nonsingular.
    TwistedHessianCurve::new_normal_form(f(8))
}

fn infinity_points(curve: &TwistedHessianCurve<F103>) -> Vec<TwistedHessianPoint<F103>> {
    let mut pts = Vec::new();

    for y in 0..103i64 {
        let p = TwistedHessianPoint::new(F103::one(), f(y), F103::zero());
        if curve.is_on_curve(&p) {
            pts.push(p);
        }
    }

    pts
}

fn all_points(curve: &TwistedHessianCurve<F103>) -> Vec<TwistedHessianPoint<F103>> {
    let mut pts = Vec::new();

    for x in 0..103i64 {
        for y in 0..103i64 {
            if curve.contains_affine(&f(x), &f(y)) {
                pts.push(TwistedHessianPoint::from_affine(f(x), f(y)));
            }
        }
    }

    pts.extend(infinity_points(curve));
    pts
}

#[test]
fn twisted_hessian_identity_is_on_curve() {
    let curve = test_curve();
    let id = TwistedHessianPoint::<F103>::identity();
    assert!(curve.is_on_curve(&id));
    assert!(id.is_identity());
}

#[test]
fn twisted_hessian_infinity_points_are_on_curve() {
    let curve = test_curve();
    let inf = infinity_points(&curve);
    assert!(!inf.is_empty());

    for p in inf {
        assert!(curve.is_on_curve(&p));
        assert!(p.is_at_infinity());
    }
}

#[test]
fn twisted_hessian_identity_is_neutral_for_all_points() {
    let curve = test_curve();
    let id = TwistedHessianPoint::<F103>::identity();

    for p in all_points(&curve) {
        assert_eq!(p.add(&id, &curve), p, "P + O != P for P={:?}", p);
        assert_eq!(id.add(&p, &curve), p, "O + P != P for P={:?}", p);
    }
}

#[test]
fn twisted_hessian_negation_works_for_all_points() {
    let curve = test_curve();
    let id = TwistedHessianPoint::<F103>::identity();

    for p in all_points(&curve) {
        let neg = p.negate(&curve);
        assert_eq!(p.add(&neg, &curve), id, "P + (-P) != O for P={:?}", p);
        assert_eq!(neg.add(&p, &curve), id, "(-P) + P != O for P={:?}", p);
    }
}

#[test]
fn twisted_hessian_double_matches_add_self() {
    let curve = test_curve();

    for p in all_points(&curve) {
        assert_eq!(
            p.double(&curve),
            p.add(&p, &curve),
            "[2]P mismatch for P={:?}",
            p
        );
    }
}

#[test]
fn twisted_hessian_addition_is_commutative() {
    let curve = test_curve();
    let pts = all_points(&curve);

    for p in &pts {
        for q in &pts {
            assert_eq!(
                p.add(q, &curve),
                q.add(p, &curve),
                "P+Q != Q+P for P={:#}, Q={:#}",
                p,
                q,
            );
        }
    }
}

#[test]
fn twisted_hessian_addition_is_associative() {
    let curve = test_curve();
    let pts = all_points(&curve);

    for p in &pts {
        for q in &pts {
            for r in &pts {
                assert_eq!(
                    p.add(q, &curve).add(r, &curve),
                    p.add(&q.add(r, &curve), &curve),
                    "(P+Q)+R != P+(Q+R) for P={:#}, Q={:#}, R={:#}",
                    p,
                    q,
                    r,
                );
            }
        }
    }
}

#[test]
fn twisted_hessian_scalar_mul_matches_repeated_addition() {
    let curve = test_curve();
    let p = TwistedHessianPoint::from_affine(f(16), f(101));
    assert!(curve.is_on_curve(&p));

    let seven_p = p.scalar_mul(&[7], &curve);
    let mut acc = TwistedHessianPoint::identity();
    for _ in 0..7 {
        acc = acc.add(&p, &curve);
    }

    assert_eq!(seven_p, acc);
}

#[test]
fn twisted_hessian_group_order_annihilates_points() {
    let curve = test_curve();
    let pts = all_points(&curve);
    let order = pts.len() as u64;

    for p in pts.iter().take(12) {
        assert_eq!(
            p.scalar_mul(&[order], &curve),
            TwistedHessianPoint::identity(),
            "[order]P != O for P={:?}",
            p,
        );
    }
}

#[test]
fn twisted_hessian_j_invariant_matches_formula() {
    let curve = test_curve();

    let d2 = <F103 as FieldOps>::square(&curve.d);
    let d3 = &curve.d * &d2;
    let num = &(&F103::from_u64(3) * &curve.d) * &(&(&F103::from_u64(8) * &curve.a) + &d3);
    let den = &d3 - &curve.a;
    let inner = &num * &den.invert().into_option().expect("nonzero denominator");
    let expected = &curve.a.invert().into_option().expect("a invertible")
        * &(&inner * &<F103 as FieldOps>::square(&inner));

    assert_eq!(curve.j_invariant(), expected);
}
