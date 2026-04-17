use crypto_bigint::{Uint, const_prime_monty_params};
use fp::fp_element::FpElement;
use fp::field_ops::FieldOps;
use ec::curve_legendre::LegendreCurve;
use ec::curve_ops::Curve;
use ec::point_legendre::LegendrePoint;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;
fn fp(n: u64) -> F19 {
    F19::from_u64(n)
}

fn curve() -> LegendreCurve<F19> {
    LegendreCurve::new(fp(3))
}

#[test]
fn legendre_curve_detects_singularity() {
    assert!(LegendreCurve::new(F19::zero()).is_singular());
    assert!(LegendreCurve::new(F19::one()).is_singular());
    assert!(!curve().is_singular());
}

#[test]
fn legendre_curve_a_invariants_are_correct() {
    let c = curve();
    let a = c.a_invariants();

    let four = F19::from(FpElement::from_u64(4));
    let three = F19::from(FpElement::from_u64(3));

    assert_eq!(a.len(), 5);
    assert_eq!(a[0], F19::zero());
    assert_eq!(a[1], -four); // -(1 + λ) = -(1 + 3) = -4
    assert_eq!(a[2], F19::zero());
    assert_eq!(a[3], F19::from(three));
    assert_eq!(a[4], F19::zero());
}

#[test]
fn legendre_curve_identity_is_on_curve() {
    let c = curve();
    assert!(c.is_on_curve(&LegendrePoint::identity()));
}

#[test]
fn legendre_curve_visible_two_torsion_is_on_curve() {
    let c = curve();
    let three = F19::from(FpElement::from_u64(3));


    assert!(c.is_on_curve(&LegendrePoint::new(F19::zero(), F19::zero())));
    assert!(c.is_on_curve(&LegendrePoint::new(F19::one(), F19::zero())));
    assert!(c.is_on_curve(&LegendrePoint::new(three, F19::zero())));
}

#[test]
fn legendre_curve_rejects_off_curve_point() {
    let c = curve();
    let two = F19::from(FpElement::from_u64(2));

    assert!(!c.is_on_curve(&LegendrePoint::new(two, two)));
}

#[test]
fn legendre_curve_j_invariant_matches_closed_formula() {
    let c = curve();
    let lambda = F19::from(FpElement::from_u64(3));
    let two_five_six = F19::from(FpElement::from_u64(256));




    let lambda_sq = <F19 as FieldOps>::square(&lambda);
    let t = lambda_sq - lambda + F19::one();
    let num = two_five_six * t * t * t;

    let den = lambda_sq * <F19 as FieldOps>::square(&(lambda - F19::one()));
    let den_inv = den.invert().into_option().unwrap();

    let expected = num * den_inv;

    assert_eq!(c.j_invariant(), expected);
}
