use crypto_bigint::{const_prime_monty_params, Uint};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;
use fp::fp_ext::{FpExt, IrreduciblePoly, TonelliShanksConstants}; // ← was missing IrreduciblePoly
use rand::RngExt;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);

type Fp19 = FpElement<Fp19Mod, 1>;

struct QuadPoly;
struct TSQuad;

impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
    fn modulus() -> [Fp19; 2] {
        [Fp19::one(), Fp19::zero()]
    }
}

impl TonelliShanksConstants<Fp19Mod, 1, 2, 1> for TSQuad {
    // Still only need 1 limb for $19^2$
    const ORDER: Uint<1> = Uint::<1>::from_u64(360);
    const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(180);
    const S: u64 = 3;
    const T: Uint<1> = Uint::<1>::from_u64(45);
    const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(22);
    const TWOSM1: Uint<1> = Uint::<1>::from_u64(4);
    fn root_of_unity() -> [FpElement<Fp19Mod, 1>; 2] {
        [Fp19::from_u64(3), Fp19::from_u64(3)]
    }
}

type F19_2 = FpExt<Fp19Mod, 1, 2, 1, QuadPoly, TSQuad>;

fn fp(n: u64) -> Fp19 {
    Fp19::from_u64(n)
}
fn el(a: u64, b: u64) -> F19_2 {
    F19_2::new([fp(a), fp(b)])
}

// -----------------------------------------------------------------------
// Structural
// -----------------------------------------------------------------------

#[test]
fn degree_is_m() {
    assert_eq!(F19_2::degree(), 2);
}

#[test]
fn characteristic_equals_base_field() {
    assert_eq!(F19_2::characteristic(), Fp19::characteristic());
    assert_eq!(F19_2::characteristic(), vec![19u64]);
}

// -----------------------------------------------------------------------
// Identity elements
// -----------------------------------------------------------------------

#[test]
fn zero_is_zero() {
    assert!(bool::from(F19_2::zero().is_zero()));
}

#[test]
fn one_is_one() {
    assert!(bool::from(F19_2::one().is_one()));
}

/*
#[test]
fn nonzero_is_not_zero() {
    assert!(bool::from(!el(1, 0).is_zero()));
    assert!(bool::from(!el(0, 1).is_zero()));
    assert!(bool::from(!el(3, 5).is_zero()));
}
*/

// -----------------------------------------------------------------------
// Addition / subtraction / negation
// -----------------------------------------------------------------------

#[test]
fn add_coefficient_wise() {
    let c = FieldOps::add(&el(3, 2), &el(1, 1));
    assert_eq!(c.coeffs[0].as_limbs()[0], 4);
    assert_eq!(c.coeffs[1].as_limbs()[0], 3);
}

#[test]
fn add_wraps_mod_p() {
    let c = FieldOps::add(&el(18, 0), &el(5, 0));
    assert_eq!(c.coeffs[0].as_limbs()[0], 4);
}

#[test]
fn sub_coefficient_wise() {
    let c = FieldOps::sub(&el(5, 4), &el(2, 1));
    assert_eq!(c.coeffs[0].as_limbs()[0], 3);
    assert_eq!(c.coeffs[1].as_limbs()[0], 3);
}

#[test]
fn negate_then_add_is_zero() {
    let a = el(7, 3);
    assert!(bool::from(FieldOps::add(&a, &a.negate()).is_zero()));
}

#[test]
fn add_zero_is_identity() {
    let a = el(5, 11);
    assert_eq!(FieldOps::add(&a, &F19_2::zero()), a);
}

#[test]
fn double_equals_add_self() {
    let a = el(3, 7);
    assert_eq!(a.double(), FieldOps::add(&a, &a));
}

// -----------------------------------------------------------------------
// Multiplication and reduction
// -----------------------------------------------------------------------

#[test]
fn mul_by_one_is_identity() {
    let a = el(5, 11);
    assert_eq!(FieldOps::mul(&a, &F19_2::one()), a);
}

#[test]
fn mul_x_squared_reduces_to_minus_one() {
    let x = el(0, 1);
    let r = FieldOps::mul(&x, &x);
    assert_eq!(r.coeffs[0].as_limbs()[0], 18);
    assert!(bool::from(r.coeffs[1].is_zero()));
}

#[test]
fn mul_concrete_values() {
    // (3+2x)(1+x) = 1 + 5x
    let r = FieldOps::mul(&el(3, 2), &el(1, 1));
    assert_eq!(r.coeffs[0].as_limbs()[0], 1);
    assert_eq!(r.coeffs[1].as_limbs()[0], 5);
}

#[test]
fn mul_commutativity() {
    let a = el(7, 3);
    let b = el(2, 11);
    assert_eq!(FieldOps::mul(&a, &b), FieldOps::mul(&b, &a));
}

#[test]
fn mul_associativity() {
    let a = el(2, 1);
    let b = el(3, 4);
    let c = el(5, 6);
    assert_eq!(
        FieldOps::mul(&FieldOps::mul(&a, &b), &c),
        FieldOps::mul(&a, &FieldOps::mul(&b, &c))
    );
}

#[test]
fn mul_distributivity() {
    let a = el(3, 2);
    let b = el(1, 5);
    let c = el(4, 7);
    assert_eq!(
        FieldOps::mul(&a, &FieldOps::add(&b, &c)),
        FieldOps::add(&FieldOps::mul(&a, &b), &FieldOps::mul(&a, &c))
    );
}

#[test]
fn square_equals_mul_self() {
    let a = el(4, 9);
    assert_eq!(a.square(), FieldOps::mul(&a, &a));
}

// -----------------------------------------------------------------------
// Inversion
// -----------------------------------------------------------------------

#[test]
fn invert_concrete_value() {
    let a = el(3, 2);
    let inv = a.invert().expect("(3+2x) is invertible in F₁₉²");
    assert_eq!(inv.coeffs[0].as_limbs()[0], 9);
    assert_eq!(inv.coeffs[1].as_limbs()[0], 13);
}

#[test]
fn one_is_not_zero_in_fp_ext() {
    let a = el(1, 0);
    assert!(!bool::from(a.is_zero()));
}

#[test]
fn invert_correctness() {
    for (a0, a1) in [(1u64, 0u64), (0, 1), (3, 2), (7, 5), (13, 18), (1, 18)] {
        let a = el(a0, a1);
        let inv = a.invert().unwrap();
        assert!(
            bool::from(FieldOps::mul(&a, &inv).is_one()),
            "inv failed for ({},{})",
            a0,
            a1
        );
    }
}

#[test]
fn invert_zero_is_none() {
    assert!(bool::from(F19_2::zero().invert().is_none()));
}

// -----------------------------------------------------------------------
// Frobenius
// -----------------------------------------------------------------------

#[test]
fn frobenius_of_base_element_is_identity() {
    let a = F19_2::from_base(fp(7));
    assert_eq!(a.frobenius(), a);
}

#[test]
fn frobenius_of_x() {
    let frob = el(0, 1).frobenius();
    assert!(bool::from(frob.coeffs[0].is_zero()));
    assert_eq!(frob.coeffs[1].as_limbs()[0], 18);
}

#[test]
fn frobenius_squared_is_identity() {
    let a = el(5, 11);
    assert_eq!(a.frobenius().frobenius(), a);
}

// -----------------------------------------------------------------------
// Norm and trace
// -----------------------------------------------------------------------

#[test]
fn norm_lies_in_base_field() {
    let n = el(3, 4).norm();
    assert_eq!(n.coeffs[0].as_limbs()[0], 6); // 9+16=25≡6 mod 19
    assert!(bool::from(n.coeffs[1].is_zero()));
}

#[test]
fn trace_lies_in_base_field() {
    let t = el(7, 5).trace();
    assert_eq!(t.coeffs[0].as_limbs()[0], 14); // 2·7=14
    assert!(bool::from(t.coeffs[1].is_zero()));
}

#[test]
fn norm_of_base_element_is_square() {
    let n = F19_2::from_base(fp(4)).norm();
    assert_eq!(n.coeffs[0].as_limbs()[0], 16); // 4²=16
}

// -----------------------------------------------------------------------
// pow
// -----------------------------------------------------------------------

#[test]
fn pow_zero_is_one() {
    assert!(bool::from(el(3, 7).pow(&[0]).is_one()));
}

#[test]
fn pow_one_is_self() {
    let a = el(3, 7);
    assert_eq!(a.pow(&[1]), a);
}

#[test]
fn pow_group_order() {
    // |F₁₉²*| = 19² − 1 = 360
    assert!(bool::from(el(3, 2).pow(&[360]).is_one()));
}

#[test]
fn pow_vartime_group_order() {
    // |F₁₉²*| = 19² − 1 = 360
    assert!(bool::from(el(3, 2).pow_vartime(&[360]).is_one()));
}

#[test]
fn pow_vartime_eq_pow() {
    // |F₁₉²*| = 19² − 1 = 360
    let x = el(3, 2).pow_vartime(&[11]);
    let y = el(3, 2).pow(&[11]);
    assert!(bool::from(x.sub(&y).is_zero()));
}

// -----------------------------------------------------------------------
// Legendre and squareroot
// -----------------------------------------------------------------------
#[test]
fn quadratic_legendre_of_qr() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut a = el(x, y);
    a = a.pow(&[2]);
    assert_eq!(a.legendre(), 1);
}

#[test]
fn ts_consts_test_quad() {
    let z = F19_2::new(TSQuad::root_of_unity());
    let mut exp = 1u64;
    for _ in 1..=TSQuad::S {
        exp *= 2;
    }
    assert_eq!(z.pow(&[exp]), F19_2::one());
}

#[test]
fn sqrt_test_quad() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut a = el(x, y);
    a = a.pow(&[2]);
    let z = a.sqrt().unwrap();
    assert_eq!(z.pow(&[2]), a);
}

#[test]
fn sqrt_test_nr_quad() {
    let a = el(5u64, 3u64); // nonresidue
    let b = el(4u64, 7u64); // nonresidue
    assert!(bool::from(a.sqrt().is_none()));
    assert!(bool::from(b.sqrt().is_none()));
}

#[test]
fn inv_and_sqrt_test_quad() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut a = el(x, y);
    a = a.pow(&[2]);
    let (myinv, mysqrt) = a.inverse_and_sqrt();
    let should_be_one = FieldOps::mul(&a, &myinv.unwrap());
    assert_eq!(mysqrt.unwrap().pow(&[2]), a);
    assert!(bool::from(should_be_one.is_one()));
}

#[test]
fn inv_and_sqrt_test_zero_quad() {
    let a = el(0u64, 0u64); // nonresidue
    let (myinv, mysqrt) = a.inverse_and_sqrt();
    assert!(bool::from(myinv.is_none()));
    assert_eq!(mysqrt.unwrap(), a);
}

#[test]
fn inv_and_sqrt_test_nr_quad() {
    let a = el(5u64, 3u64); // nonresidue
    let (myinv, mysqrt) = a.inverse_and_sqrt();
    assert!(bool::from(a.mul(&myinv.unwrap()).is_one()));
    assert!(bool::from(mysqrt.is_none()));

    let b = el(4u64, 7u64); // nonresidue
    let (myinv, mysqrt) = b.inverse_and_sqrt();
    assert!(bool::from(b.mul(&myinv.unwrap()).is_one()));
    assert!(bool::from(mysqrt.is_none()));
}

#[test]
fn invme_sqrtother_test_quad() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let a = el(x, y);
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut b = el(x, y);
    b = b.pow(&[2]);
    let (myinv, mysqrt) = a.invertme_sqrtother(&b);
    let should_be_one = FieldOps::mul(&a, &myinv.unwrap());
    assert_eq!(mysqrt.unwrap().pow(&[2]), b);
    assert!(bool::from(should_be_one.is_one()));
}

#[test]
fn invertme_sqrtother_zero_quad() {
    let mut rng = rand::rng();
    let a = el(0u64, 0u64);
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut b = el(x, y);
    b = b.pow(&[2]);
    let (myinv, mysqrt) = a.invertme_sqrtother(&b);
    assert!(bool::from(myinv.is_none()));
    assert!(bool::from(mysqrt.is_none()));
}

#[test]
fn invertme_sqrtother_nr_quad() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let a = el(x, y);
    let b = el(5u64, 3u64); // nonresidue
    let (myinv, mysqrt) = a.invertme_sqrtother(&b);
    assert!(bool::from(a.mul(&myinv.unwrap()).is_one()));
    assert!(bool::from(mysqrt.is_none()));
}

#[test]
fn sqrtratio_test_quad() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut a = el(x, y);
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let mut b = el(x, y);
    a = a.pow(&[2]);
    b = b.pow(&[2]);
    let mysqrt = a.sqrt_ratio(&b).unwrap();
    assert_eq!(mysqrt.pow(&[2]), a.mul(&b.invert().unwrap()));
}

#[test]
fn sqrtratio_nr_quad() {
    let a = el(1, 0);
    let b = el(5u64, 3u64); // nonresidue
    let mysqrt = a.sqrt_ratio(&b);
    assert!(bool::from(mysqrt.is_none()));
}

#[test]
fn sqrtratio_zero_quad() {
    let a = el(1, 0);
    let b = el(0, 0); // nonresidue
    let mysqrt = a.sqrt_ratio(&b);
    assert!(bool::from(mysqrt.is_none()));
}

// -----------------------------------------------------------------------
// Cubic extension  F₁₉³  with  f(x) = x³ − 2  (≡ x³ + 17 mod 19)
// -----------------------------------------------------------------------

struct CubicPoly;
struct TSCubic;

impl IrreduciblePoly<Fp19Mod, 1, 3> for CubicPoly {
    fn modulus() -> [Fp19; 3] {
        [fp(17), fp(0), fp(0)] // [c₀=17, c₁=0, c₂=0]  →  x³ + 17
    }
}

impl TonelliShanksConstants<Fp19Mod, 1, 3, 1> for TSCubic {
    // Still only need 1 limb for 19^3
    const ORDER: Uint<1> = Uint::<1>::from_u64(6858);
    const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(3429);
    const S: u64 = 1;
    const T: Uint<1> = Uint::<1>::from_u64(3429);
    const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(1714);
    const TWOSM1: Uint<1> = Uint::<1>::from_u64(1);
    fn root_of_unity() -> [FpElement<Fp19Mod, 1>; 3] {
        [fp(1), fp(0), fp(0)]
    }
}

type F19_3 = FpExt<Fp19Mod, 1, 3, 1, CubicPoly, TSCubic>;

fn el3(a: u64, b: u64, c: u64) -> F19_3 {
    F19_3::new([fp(a), fp(b), fp(c)])
}

#[test]
fn cubic_degree_is_3() {
    assert_eq!(F19_3::degree(), 3);
}
#[test]
fn cubic_zero_one() {
    assert!(bool::from(F19_3::zero().is_zero()));
    assert!(bool::from(F19_3::one().is_one()));
}

#[test]
fn cubic_x_cubed_reduces() {
    // x³ ≡ 2 mod (x³−2)
    let x = el3(0, 1, 0);
    let x3 = FieldOps::mul(&FieldOps::mul(&x, &x), &x);
    assert_eq!(x3.coeffs[0].as_limbs()[0], 2);
    assert!(bool::from(x3.coeffs[1].is_zero()));
    assert!(bool::from(x3.coeffs[2].is_zero()));
}

#[test]
fn cubic_mul_commutativity() {
    let a = el3(1, 2, 3);
    let b = el3(4, 5, 6);
    assert_eq!(FieldOps::mul(&a, &b), FieldOps::mul(&b, &a));
}

#[test]
fn cubic_invert_correctness() {
    let a = el3(1, 2, 3);
    assert!(bool::from(
        FieldOps::mul(&a, &a.invert().expect("invertible")).is_one()
    ));
}

#[test]
fn cubic_pow_order() {
    // |F₁₉³*| = 19³ − 1 = 6858
    assert!(bool::from(el3(2, 1, 5).pow(&[6858]).is_one()));
}

#[test]
fn ts_consts_test_cubic() {
    let z = F19_3::new(TSCubic::root_of_unity());
    assert_eq!(z.pow(&[TSCubic::S]), F19_3::one());
}

#[test]
fn sqrt_test_cubic() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let z: u64 = rng.random();
    let mut a = el3(x, y, z);
    a = a.pow(&[2]);
    let z = a.sqrt().unwrap();
    assert_eq!(z.pow(&[2]), a);
}

//////////////
#[test]
fn sqrt_test_nr_cubic() {
    let a = el3(1u64, 2u64, 3u64); // nonresidue
    let b = el3(4u64, 4u64, 3u64); // nonresidue
    assert!(bool::from(a.sqrt().is_none()));
    assert!(bool::from(b.sqrt().is_none()));
}

#[test]
fn inv_and_sqrt_test_cubic() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let z: u64 = rng.random();
    let mut a = el3(x, y, z);
    a = a.pow(&[2]);
    let (myinv, mysqrt) = a.inverse_and_sqrt();
    let should_be_one = FieldOps::mul(&a, &myinv.unwrap());
    assert_eq!(mysqrt.unwrap().pow(&[2]), a);
    assert!(bool::from(should_be_one.is_one()));
}

#[test]
fn inv_and_sqrt_test_zero_cubic() {
    let a = el3(0u64, 0u64, 0u64);
    let (myinv, mysqrt) = a.inverse_and_sqrt();
    assert!(bool::from(myinv.is_none()));
    assert_eq!(mysqrt.unwrap(), a);
}

#[test]
fn inv_and_sqrt_test_nr_cubic() {
    let a = el3(1u64, 2u64, 3u64); // nonresidue
    let (myinv, mysqrt) = a.inverse_and_sqrt();
    assert!(bool::from(a.mul(&myinv.unwrap()).is_one()));
    assert!(bool::from(mysqrt.is_none()));
}

#[test]
fn invme_sqrtother_test_cubic() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let z: u64 = rng.random();
    let a = el3(x, y, z);
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let z: u64 = rng.random();
    let mut b = el3(x, y, z);
    b = b.pow(&[2]);
    let (myinv, mysqrt) = a.invertme_sqrtother(&b);
    let should_be_one = FieldOps::mul(&a, &myinv.unwrap());
    assert_eq!(mysqrt.unwrap().pow(&[2]), b);
    assert!(bool::from(should_be_one.is_one()));
}

#[test]
fn invertme_sqrtother_nr_cubic() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let z: u64 = rng.random();
    let a = el3(x, y, z);
    let b = el3(1u64, 2u64, 3u64); // nonresidue
    let (myinv, mysqrt) = a.invertme_sqrtother(&b);
    assert!(bool::from(a.mul(&myinv.unwrap()).is_one()));
    assert!(bool::from(mysqrt.is_none()));
}

#[test]
fn sqrtratio_test_cubic() {
    let mut rng = rand::rng();
    let x: u64 = rng.random();
    let y: u64 = rng.random();
    let z: u64 = rng.random();
    let mut a = el3(x, y, z);
    let x: u64 = rng.random();
    let z: u64 = rng.random();
    let mut b = el3(x, y, z);
    a = a.pow(&[2]);
    b = b.pow(&[2]);
    let mysqrt = a.sqrt_ratio(&b).unwrap();
    assert_eq!(mysqrt.pow(&[2]), a.mul(&b.invert().unwrap()));
}

#[test]
fn sqrtratio_nr_cubic() {
    let a = el3(1, 0, 0);
    let b = el3(1u64, 2u64, 3u64); // nonresidue
    let mysqrt = a.sqrt_ratio(&b);
    assert!(bool::from(mysqrt.is_none()));
}
