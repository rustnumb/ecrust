use crypto_bigint::{const_prime_monty_params, Uint};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;
use fp::fp_ext::{FpExt, IrreduciblePoly, MultiplicativeGroupOrder}; // ← was missing IrreduciblePoly

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);

type Fp19 = FpElement<Fp19Mod, 1>;

struct QuadPoly;
struct OrderQuad;

impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
    fn modulus() -> [Fp19; 2] {
        [Fp19::one(), Fp19::zero()]
    }
}

impl MultiplicativeGroupOrder<Fp19Mod, 1, 2, 1> for OrderQuad {
    // Still only need 1 limb for 19^2
    type Order = Uint<1>;

    fn order() -> Self::Order {
        const ORDER: Uint<1> = Uint::<1>::from_u64(360);
        ORDER
    }

    fn half_order() -> Self::Order {
        const ORDER: Uint<1> = Uint::<1>::from_u64(360) >> 1;
        ORDER
    }
}

type F19_2 = FpExt<Fp19Mod, 1, 2, QuadPoly, OrderQuad>;

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
    assert!(F19_2::zero().is_zero());
}

#[test]
fn one_is_one() {
    assert!(F19_2::one().is_one());
}

#[test]
fn nonzero_is_not_zero() {
    assert!(!el(1, 0).is_zero());
    assert!(!el(0, 1).is_zero());
    assert!(!el(3, 5).is_zero());
}

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
    assert!(FieldOps::add(&a, &a.negate()).is_zero());
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
    assert!(r.coeffs[1].is_zero());
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
fn invert_correctness() {
    for (a0, a1) in [(1u64, 0u64), (0, 1), (3, 2), (7, 5), (13, 18), (1, 18)] {
        let a = el(a0, a1);
        let inv = a
            .invert()
            .unwrap_or_else(|| panic!("({},{}) not invertible", a0, a1));
        assert!(
            FieldOps::mul(&a, &inv).is_one(),
            "inv failed for ({},{})",
            a0,
            a1
        );
    }
}

#[test]
fn invert_zero_is_none() {
    assert!(F19_2::zero().invert().is_none());
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
    assert!(frob.coeffs[0].is_zero());
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
    assert!(n.coeffs[1].is_zero());
}

#[test]
fn trace_lies_in_base_field() {
    let t = el(7, 5).trace();
    assert_eq!(t.coeffs[0].as_limbs()[0], 14); // 2·7=14
    assert!(t.coeffs[1].is_zero());
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
    assert!(el(3, 7).pow(&[0]).is_one());
}
#[test]
fn pow_one_is_self() {
    let a = el(3, 7);
    assert_eq!(a.pow(&[1]), a);
}

#[test]
fn pow_group_order() {
    // |F₁₉²*| = 19² − 1 = 360
    assert!(el(3, 2).pow(&[360]).is_one());
}

// -----------------------------------------------------------------------
// Cubic extension  F₁₉³  with  f(x) = x³ − 2  (≡ x³ + 17 mod 19)
// -----------------------------------------------------------------------

struct CubicPoly;
struct OrderCubic;

impl IrreduciblePoly<Fp19Mod, 1, 3> for CubicPoly {
    fn modulus() -> [Fp19; 3] {
        [fp(17), fp(0), fp(0)] // [c₀=17, c₁=0, c₂=0]  →  x³ + 17
    }
}

impl MultiplicativeGroupOrder<Fp19Mod, 1, 2, 1> for OrderCubic {
    // Still only need 1 limb for 19^3
    type Order = Uint<1>;
    fn order() -> Self::Order {
        const ORDER: Uint<1> = Uint::<1>::from_u64(6858);
        ORDER
    }
    fn half_order() -> Self::Order {
        const ORDER: Uint<1> = Uint::<1>::from_u64(6858) >> 1;
        ORDER
    }
}

type F19_3 = FpExt<Fp19Mod, 1, 3, CubicPoly, OrderCubic>;

fn el3(a: u64, b: u64, c: u64) -> F19_3 {
    F19_3::new([fp(a), fp(b), fp(c)])
}

#[test]
fn cubic_degree_is_3() {
    assert_eq!(F19_3::degree(), 3);
}
#[test]
fn cubic_zero_one() {
    assert!(F19_3::zero().is_zero());
    assert!(F19_3::one().is_one());
}

#[test]
fn cubic_x_cubed_reduces() {
    // x³ ≡ 2 mod (x³−2)
    let x = el3(0, 1, 0);
    let x3 = FieldOps::mul(&FieldOps::mul(&x, &x), &x);
    assert_eq!(x3.coeffs[0].as_limbs()[0], 2);
    assert!(x3.coeffs[1].is_zero());
    assert!(x3.coeffs[2].is_zero());
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
    assert!(FieldOps::mul(&a, &a.invert().expect("invertible")).is_one());
}

#[test]
fn cubic_pow_order() {
    // |F₁₉³*| = 19³ − 1 = 6858
    assert!(el3(2, 1, 5).pow(&[6858]).is_one());
}
