use crypto_bigint::{const_prime_monty_params, Uint};
use ec::curve_legendre::LegendreCurve;
use ec::curve_ops::Curve;
use ec::point_legendre::LegendrePoint;
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

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
    assert_eq!(a[1], -&four); // -(1 + λ) = -(1 + 3) = -4
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
    println!("lambda = {}", lambda.as_uint().to_words()[0]);
    let two_five_six = F19::from(FpElement::from_u64(256));

    let lambda_sq = <F19 as FieldOps>::square(&lambda);
    let t = &(&lambda_sq - &lambda) + &F19::one();
    println!("t = {}", t.as_uint().to_words()[0]);
    let num = &(&two_five_six * &t) * &(&t * &t);
    println!("num = {}", num.as_uint().to_words()[0]);

    let den = &lambda_sq * &<F19 as FieldOps>::square(&(&lambda - &F19::one()));
    println!("den = {}", den.as_uint().to_words()[0]);
    let den_inv = den.invert().into_option().unwrap();

    let expected = &num * &den_inv;
    println!("expected = {}", expected.as_uint().to_words()[0]);

    assert_eq!(c.j_invariant(), expected);
}

#[test]
fn legendre_to_general_weierstrass_matches_a_invariants() {
    let c = curve(); // e.g. λ = 3 over F_19
    let w = c.to_weierstrass();

    assert_eq!(w.a1, F19::zero());
    assert_eq!(w.a2, -&(&F19::one() + &c.lambda));
    assert_eq!(w.a3, F19::zero());
    assert_eq!(w.a4, c.lambda);
    assert_eq!(w.a6, F19::zero());
}

#[test]
fn legendre_to_general_weierstrass_preserves_j() {
    let c = curve();
    let w = c.to_weierstrass();

    assert_eq!(c.j_invariant(), w.j_invariant());
}

#[test]
fn legendre_to_short_weierstrass_preserves_j() {
    let c = curve();
    let w = c.to_short_weierstrass();

    assert_eq!(c.j_invariant(), w.j_invariant());
}

#[test]
fn legendre_short_weierstrass_coordinate_shift_works() {
    let c = curve();
    let jc = c.j_invariant();
    let w = c.to_short_weierstrass();
    let jw = w.j_invariant();
    assert_eq!(jc, jw);

    let p = LegendrePoint::new(fp(2), fp(6));
    assert!(c.is_on_curve(&p));

    let three = fp(3);
    let three_inv = three.invert().into_option().unwrap();
    let shift = &(&F19::one() + &c.lambda) * &three_inv;

    let xw = &p.x - &shift;
    let yw = p.y;

    assert!(w.contains(&xw, &yw));
}
