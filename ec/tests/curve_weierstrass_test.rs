// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------


use crypto_bigint::{Uint, const_prime_monty_params};
use fp::fp_element::FpElement;
use fp::field_ops::FieldOps;
use ec::curve_weierstrass::WeierstrassCurve;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(n: u64) -> F19 {
    F19::from_u64(n)
}

#[test]
fn short_weierstrass_has_zero_a1_a2_a3() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    assert!(bool::from(c.a1.is_zero()));
    assert!(bool::from(c.a2.is_zero()));
    assert!(bool::from(c.a3.is_zero()));
    assert_eq!(c.a4, fp(2));
    assert_eq!(c.a6, fp(3));
}

#[test]
fn contains_known_point() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    let rhs = fp(6);
    let y = rhs.sqrt().into_option().expect("6 must be a QR mod 19");
    assert!(c.contains(&fp(1), &y));
}

#[test]
fn rejects_non_point() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    assert!(!c.contains(&fp(0), &fp(0)));
}

#[test]
fn non_singular_curve() {
    let c = WeierstrassCurve::new_short(fp(2), fp(3));
    assert!(!bool::from(c.discriminant().is_zero()));
}
