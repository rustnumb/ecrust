//! Integration tests for generalized Hessian curves.

use crypto_bigint::{Uint, const_prime_monty_params};

use ec::curve_hessian::HessianCurve;
use ec::curve_ops::Curve;
use ec::point_hessian::HessianPoint;
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp103Mod, Uint<1>, "0000000000000067", 2);
type F103 = FpElement<Fp103Mod, 1>;

fn f(n: i64) -> F103 {
    F103::from_u64(n.rem_euclid(103) as u64)
}

fn cube_roots_of_unity() -> [F103; 3] {
    [f(1), f(46), f(56)]
}

fn complete_curve() -> HessianCurve<F103> {
    // Over F_103, 2 is not a cube, so H_{2,1} is complete by Theorem 3.
    HessianCurve::new(f(2), f(1))
}

fn all_points(curve: &HessianCurve<F103>) -> Vec<HessianPoint<F103>> {
    let mut pts = Vec::new();

    for x in 0..103i64 {
        for y in 0..103i64 {
            if curve.contains_affine(&f(x), &f(y)) {
                pts.push(HessianPoint::from_affine(f(x), f(y)));
            }
        }
    }

    for omega in cube_roots_of_unity() {
        pts.push(HessianPoint::new(F103::one(), -omega, F103::zero()));
    }

    pts
}

#[test]
fn hessian_identity_is_on_curve() {
    let curve = complete_curve();
    let id = HessianPoint::<F103>::identity();
    assert!(curve.is_on_curve(&id));
    assert!(id.is_identity());
}

#[test]
fn hessian_infinity_points_are_on_curve() {
    let curve = complete_curve();
    for omega in cube_roots_of_unity() {
        let p = HessianPoint::new(F103::one(), -omega, F103::zero());
        assert!(curve.is_on_curve(&p));
        assert!(p.is_at_infinity());
    }
}

#[test]
fn hessian_identity_is_neutral_for_all_points() {
    let curve = complete_curve();
    let id = HessianPoint::<F103>::identity();

    for p in all_points(&curve) {
        assert_eq!(p.add(&id, &curve), p, "P + O != P for P={:#}", p);
        assert_eq!(id.add(&p, &curve), p, "O + P != P for P={:#}", p);
    }
}

#[test]
fn hessian_negation_works_for_all_points() {
    let curve = complete_curve();
    let id = HessianPoint::<F103>::identity();

    for p in all_points(&curve) {
        let neg = p.negate(&curve);
        assert_eq!(p.add(&neg, &curve), id, "P + (-P) != O for P={:#}", p);
        assert_eq!(neg.add(&p, &curve), id, "(-P) + P != O for P={:#}", p);
    }
}

#[test]
fn hessian_addition_is_commutative_on_complete_curve() {
    let curve = complete_curve();
    let pts = all_points(&curve);

    for p in &pts {
        for q in &pts {
            assert_eq!(
                p.add(q, &curve),
                q.add(p, &curve),
                "P+Q != Q+P for P={:#}, Q={:#}",
                p,
                q
            );
        }
    }
}

#[test]
fn hessian_addition_is_associative_on_complete_curve() {
    let curve = complete_curve();
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
fn hessian_scalar_mul_matches_repeated_addition() {
    let curve = complete_curve();
    let p = HessianPoint::from_affine(f(16), f(11));
    assert!(curve.is_on_curve(&p));

    let seven_p = p.scalar_mul(&[7], &curve);
    let mut acc = HessianPoint::identity();
    for _ in 0..7 {
        acc = acc.add(&p, &curve);
    }

    assert_eq!(seven_p, acc);
}

#[test]
fn hessian_group_order_annihilates_points() {
    let curve = complete_curve();
    let pts = all_points(&curve);
    let order = pts.len() as u64;

    for p in pts.iter().take(12) {
        assert_eq!(p.scalar_mul(&[order], &curve), HessianPoint::identity());
    }
}

#[test]
fn hessian_efd_constructor_matches_paper_parameterization() {
    let via_efd = HessianCurve::new_efd(f(5));
    let via_paper = HessianCurve::new_hessian(f(15));
    assert_eq!(via_efd, via_paper);
}

#[test]
fn hessian_weierstrass_birational_roundtrip() {
    let curve = HessianCurve::new_hessian(f(15));
    let zeta = F103::one();

    let p = HessianPoint::from_affine(f(3), f(24));
    assert!(curve.is_on_curve(&p));

    let wc = curve.to_weierstrass_curve_with_zeta(zeta).unwrap();

    let wp = curve
        .map_point_to_weierstrass_with_zeta(&p, zeta)
        .expect("forward birational map");
    assert!(wc.is_on_curve(&wp));

    let back = curve
        .map_point_from_weierstrass_with_zeta(&wp, zeta)
        .expect("inverse birational map");
    assert!(curve.is_on_curve(&back));

    assert_eq!(p, back);
}