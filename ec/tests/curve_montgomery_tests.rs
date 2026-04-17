use crypto_bigint::{Uint, const_prime_monty_params};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

use ec::curve_montgomery::MontgomeryCurve;
use ec::curve_ops::Curve;
use ec::point_montgomery::KummerPoint;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(n: u64) -> F19 {
    F19::from_u64(n)
}

/// Brute-force: all affine x-coordinates for which
///
///     y² = (x³ + A x² + x)/B
///
/// has a solution over F_19.
fn all_x_coords_montgomery_19(a: F19, b: F19) -> Vec<F19> {
    let mut xs = Vec::new();

    for xint in 0..19u64 {
        let x = fp(xint);
        let x2 = &x * &x;
        let ax = &a * &x;
        let rhs_num = &x * &(&x2 + &(&ax + &F19::one()));
        let rhs = &rhs_num * &F19::invert(&b).expect("b should be invertible mod 19");

        let mut has_y = false;
        for y in 0..19u64 {
            if fp(y * y) == rhs {
                has_y = true;
                break;
            }
        }

        if has_y {
            xs.push(x);
        }
    }

    xs
}

#[test]
fn montgomery_new_stores_coefficients() {
    let c = MontgomeryCurve::new(fp(3), fp(5));
    assert_eq!(c.a, fp(3));
    assert_eq!(c.b, fp(5));
}

#[test]
fn montgomery_a_invariants_are_ab() {
    let c = MontgomeryCurve::new(fp(3), fp(5));
    let inv = c.a_invariants();
    assert_eq!(inv, [fp(3), fp(5)]);
}

#[test]
fn non_singular_curve() {
    assert!(MontgomeryCurve::<F19>::is_smooth(&fp(3), &fp(5)));
    assert!(MontgomeryCurve::<F19>::is_smooth(&fp(7), &fp(1)));
}

#[test]
fn montgomery_curve_rejects_b_zero() {
    assert!(!MontgomeryCurve::<F19>::is_smooth(&fp(3), &fp(0)));
}

#[test]
fn montgomery_curve_rejects_a_equal_plus_minus_two() {
    let two = fp(2);
    let minus_two = -&two;

    assert!(!MontgomeryCurve::<F19>::is_smooth(&two, &fp(1)));
    assert!(!MontgomeryCurve::<F19>::is_smooth(&minus_two, &fp(1)));
}

#[test]
fn montgomery_a24_matches_formula() {
    let c = MontgomeryCurve::new(fp(3), fp(5));

    let four_inv = fp(4).invert().unwrap();
    let expected = &(&fp(3) + &fp(2)) * &four_inv;

    assert_eq!(c.a24(), expected);
}

#[test]
fn montgomery_known_x_coordinates_are_on_curve() {
    let c = MontgomeryCurve::new(fp(3), fp(5));

    for x in all_x_coords_montgomery_19(fp(3), fp(5)) {
        let p = KummerPoint::from_x(x);
        assert!(c.is_on_curve(&p), "x = {:?} should lie on the Montgomery curve", x);
    }
}

#[test]
fn montgomery_rejects_some_non_points() {
    let c = MontgomeryCurve::new(fp(3), fp(5));

    let all_x = all_x_coords_montgomery_19(fp(3), fp(5));
    for x in 0..19u64 {
        println!("x = {}", x);
        let p = KummerPoint::from_x(fp(x));
        let expected = all_x.contains(&fp(x));
        assert_eq!(
            c.is_on_curve(&p),
            expected,
            "unexpected on-curve result for x = {}",
            x
        );
    }
}

#[test]
fn montgomery_j_invariant_matches_formula() {
    let c = MontgomeryCurve::new(fp(3), fp(5));
    let a = fp(3);

    let asq = a.square();
    let num = &fp(256) * &(&(&asq - &fp(3)) * &(&asq - &fp(3)).square());
    let den_inv = (&asq - &fp(4)).invert().unwrap();
    let expected = &num * &den_inv;

    assert_eq!(c.j_invariant(), expected);
}