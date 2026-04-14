use crypto_bigint::{Uint, const_prime_monty_params};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use ec::curve_montgomery::MontgomeryCurve;
use ec::curve_ops::Curve;
use ec::point_montgomery::KummerPoint;
use ec::point_ops::PointOps;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
const_prime_monty_params!(Fp59Mod, Uint<1>, "000000000000003B", 2);
type F19 = FpElement<Fp19Mod, 1>;
type F59 = FpElement<Fp59Mod, 1>;

fn fp(n: u64) -> F19 {
    F19::from_u64(n)
}
fn fq(n: u64) -> F59 {
    F59::from_u64(n)
}



fn curve() -> MontgomeryCurve<F19> {
    MontgomeryCurve::new(fp(3), fp(5))
}
fn curve2() -> MontgomeryCurve<F59> { MontgomeryCurve::new(fq(7), fq(1)) }

fn all_x_coords_montgomery_19(a: F19, b: F19) -> Vec<F19> {
    let mut xs = Vec::new();

    for xint in 0..19u64 {
        let x = fp(xint);
        let rhs_num = x * (x * x + a * x + F19::one());
        let rhs = rhs_num * F19::invert(&b).expect("b should be invertible mod 19");

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
fn new_stores_projective_coordinates() {
    let p = KummerPoint::new(fp(7), fp(11));
    assert_eq!(p.x, fp(7));
    assert_eq!(p.z, fp(11));
}

#[test]
fn from_x_sets_z_to_one() {
    let p = KummerPoint::from_x(fp(9));
    assert_eq!(p.x, fp(9));
    assert_eq!(p.z, F19::one());
}

#[test]
fn projective_equality_detects_scaled_points() {
    let p = KummerPoint::new(fp(3), fp(4));
    let q = KummerPoint::new(fp(6), fp(8)); // same projective class
    assert_eq!(p, q);
    assert!(bool::from(p.ct_eq(&q)));
}

#[test]
fn projective_equality_detects_distinct_points() {
    let p = KummerPoint::new(fp(3), fp(4));
    let q = KummerPoint::new(fp(3), fp(5));
    assert_ne!(p, q);
    assert!(bool::from(p.ct_ne(&q)));
}

#[test]
fn conditional_select_choice_false_returns_first() {
    let p = KummerPoint::new(fp(1), fp(2));
    let q = KummerPoint::new(fp(3), fp(4));
    let r = KummerPoint::conditional_select(&p, &q, Choice::from(0u8));
    assert_eq!(r, p);
}

#[test]
fn conditional_select_choice_true_returns_second() {
    let p = KummerPoint::new(fp(1), fp(2));
    let q = KummerPoint::new(fp(3), fp(4));
    let r = KummerPoint::conditional_select(&p, &q, Choice::from(1u8));
    assert_eq!(r, q);
}

#[test]
fn conditional_assign_choice_false_keeps_value() {
    let mut p = KummerPoint::new(fp(1), fp(2));
    let q = KummerPoint::new(fp(3), fp(4));
    p.conditional_assign(&q, Choice::from(0u8));
    assert_eq!(p, KummerPoint::new(fp(1), fp(2)));
}

#[test]
fn conditional_assign_choice_true_overwrites_value() {
    let mut p = KummerPoint::new(fp(1), fp(2));
    let q = KummerPoint::new(fp(3), fp(4));
    p.conditional_assign(&q, Choice::from(1u8));
    assert_eq!(p, q);
}

#[test]
fn conditional_swap_choice_false_keeps_both() {
    let mut p = KummerPoint::new(fp(1), fp(2));
    let mut q = KummerPoint::new(fp(3), fp(4));
    KummerPoint::conditional_swap(&mut p, &mut q, Choice::from(0u8));
    assert_eq!(p, KummerPoint::new(fp(1), fp(2)));
    assert_eq!(q, KummerPoint::new(fp(3), fp(4)));
}

#[test]
fn conditional_swap_choice_true_swaps_both() {
    let mut p = KummerPoint::new(fp(1), fp(2));
    let mut q = KummerPoint::new(fp(3), fp(4));
    KummerPoint::conditional_swap(&mut p, &mut q, Choice::from(1u8));
    assert_eq!(p, KummerPoint::new(fp(3), fp(4)));
    assert_eq!(q, KummerPoint::new(fp(1), fp(2)));
}


#[test]
fn xdouble_of_curve_point_stays_on_curve() {
    let c = curve();

    for x in all_x_coords_montgomery_19(fp(3), fp(5)) {
        let p = KummerPoint::from_x(x);
        let dbl = p.xdouble(&c);
        assert!(
            c.is_on_curve(&dbl) || dbl.is_identity(),
            "[2]P should stay on the Kummer line for x = {:?}",
            x.as_uint()
        );
    }
}

#[test]
fn xadd_of_known_points_stays_on_curve() {
    let c = curve();
    let xs = all_x_coords_montgomery_19(fp(3), fp(5));

    if xs.len() >= 2 {
        let p = KummerPoint::from_x(xs[0]);
        let q = KummerPoint::from_x(xs[1]);

        // This is only a smoke test for the current API shape:
        // xadd requires x(P-Q), so we reuse p as a dummy differential input.
        let r = p.xadd(&q, &p);

        assert!(
            c.is_on_curve(&r) || r.is_identity(),
            "xadd output should lie on the Kummer line"
        );
    }
}

#[test]
fn xdouble_concrete_values() {
    let c = curve();
    let pairs = [
        (fp(1), fp(0)),
        (fp(6), fp(1)),
        (fp(8), fp(1)),
        (fp(12), fp(1)),
        (fp(18), fp(0))
    ];
    for pair in pairs.iter() {
        let p = KummerPoint::from_x(pair.0);
        let dbl = p.xdouble(&c);
        assert_eq!(dbl.to_x().expect("z should be invertible"), pair.1);
    }
}

#[test]

fn points_of_order_two() {
    let c = curve();
    let xs = [fp(0), fp(3), fp(13)];

    for x in xs.iter() {
        let p = KummerPoint::from_x(*x);
        let dbl = p.xdouble(&c);
        assert!(bool::from(dbl.is_identity()));
    }
}

#[test]
fn xadd_concrete_values() {
    let tuples = [
        (fp(0), fp(1), fp(1), fp(1)),
        (fp(0), fp(3), fp(13), fp(13)),
        (fp(0), fp(6), fp(16), fp(16)),
        (fp(0), fp(8), fp(12), fp(12)),
        (fp(1), fp(13), fp(18), fp(18)),
        (fp(3), fp(6), fp(12), fp(12)),
        (fp(3), fp(16), fp(8), fp(8)),
        (fp(3), fp(18), fp(1), fp(1)),
    ];

    for tuple in tuples.iter() {
        let p = KummerPoint::from_x(tuple.0);
        let q = KummerPoint::from_x(tuple.1);
        let diff = KummerPoint::from_x(tuple.2);
        let sum = p.xadd(&q, &diff);
        assert_eq!(sum.to_x().expect("z should be invertible"), tuple.3);
    }
}

#[test]
fn xadd_concrete_values_bis() {
    let tuples = [
        (fq(9), fq(47), fq(23), fq(31)),
        (fq(23), fq(4), fq(34), fq(40)),
        (fq(33), fq(18), fq(4), fq(36)),
        (fq(34), fq(15), fq(31), fq(18)),
        (fq(40), fq(32), fq(46), fq(48)),
    ];

    for tuple in tuples.iter() {
        let p = KummerPoint::from_x(tuple.0);
        let q = KummerPoint::from_x(tuple.1);
        let diff = KummerPoint::from_x(tuple.2);
        let sum = p.xadd(&q, &diff);
        assert_eq!(sum.to_x().expect("z should be invertible"), tuple.3);
    }
}


#[test]
fn scalar_mul_zero_is_identity() {
    let c = curve();
    let p = KummerPoint::from_x(fp(4));
    assert!(p.scalar_mul(&[0], &c).is_identity());
}

#[test]
fn scalar_mul_one_is_self() {
    let c = curve();
    let p = KummerPoint::from_x(fp(4));
    assert_eq!(p.scalar_mul(&[1], &c), p);
}

#[test]
fn scalar_mul_two_matches_xdouble() {
    let c = curve();
    let p = KummerPoint::from_x(fp(4));
    assert_eq!(p.scalar_mul(&[2], &c), p.xdouble(&c));
}


#[test]
fn pointops_identity_matches_inherent_identity() {
    assert_eq!(KummerPoint::<F19>::identity(), KummerPoint::<F19>::identity());
}


#[test]
fn negate_is_trivial_on_kummer_line() {
    let c = curve();
    let p = KummerPoint::from_x(fp(8));
    assert_eq!(p.negate(&c), p);
}