//! Integration tests for the `fp` crate — FpElement over a concrete prime.
//!
//! Integration tests in `tests/` are compiled as a *separate crate* that
//! links against `fp`.  The rules differ from unit tests inside `src/`:
//!   - NO `#[cfg(test)]` wrapper  (the whole file is already test-only)
//!   - NO `mod tests { }` block   (no module nesting needed)
//!   - Trait methods require the trait to be in scope with `use`

use crypto_bigint::{const_prime_monty_params, Uint};

use fp::field_ops::FieldOps;        // brings zero/one/is_zero/invert/… into scope
use fp::fp_element::FpElement;      // the concrete type

// ---------------------------------------------------------------------------
// Define the modulus p = 19 at compile time.
//
// const_prime_monty_params!(TypeName, UintType, hex_modulus, label)
// ---------------------------------------------------------------------------
const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2);

type F19 = FpElement<Fp19Modulus, 1>;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn zero_is_zero() {
    assert!(F19::zero().is_zero());
}

#[test]
fn one_is_one() {
    assert!(F19::one().is_one());
}

#[test]
fn degree_of_base_field_is_one() {
    assert_eq!(F19::degree(), 1);
}

#[test]
fn add_mod_p() {
    // 17 + 5 = 22 ≡ 3  (mod 19)
    let a = F19::from_u64(17);
    let b = F19::from_u64(5);
    assert_eq!((a + b).as_limbs()[0], 3);
}

#[test]
fn sub_mod_p() {
    // 3 − 7 = −4 ≡ 15  (mod 19)
    let a = F19::from_u64(3);
    let b = F19::from_u64(7);
    assert_eq!((a - b).as_limbs()[0], 15);
}

#[test]
fn mul_mod_p() {
    // 7 × 8 = 56 ≡ 18  (mod 19)
    let a = F19::from_u64(7);
    let b = F19::from_u64(8);
    assert_eq!((a * b).as_limbs()[0], 18);
}

#[test]
fn neg_mod_p() {
    // −3 ≡ 16  (mod 19)
    let a = F19::from_u64(3);
    assert_eq!((-a).as_limbs()[0], 16);
}

#[test]
fn square() {
    // 4² = 16  (mod 19)
    let a = F19::from_u64(4);
    assert_eq!(a.square().as_limbs()[0], 16);
}

#[test]
fn double() {
    // 2 × 9 = 18  (mod 19)
    let a = F19::from_u64(9);
    assert_eq!(a.double().as_limbs()[0], 18);
}

#[test]
fn inv_works() {
    // 7 × 7⁻¹ ≡ 1  (mod 19)
    let a   = F19::from_u64(7);
    let inv = a.invert().unwrap();
    assert_eq!((a * inv).as_limbs()[0], 1);
}

#[test]
fn inv_zero_is_none() {
    assert!(F19::zero().invert().is_none());
}

#[test]
fn characteristic_is_p() {
    // characteristic should be [19] as little-endian u64 words
    assert_eq!(F19::characteristic(), vec![19u64]);
}

#[test]
fn pow_works() {
    // 2^10 = 1024 ≡ 1024 − 53×19 = 1024 − 1007 = 17  (mod 19)
    let a = F19::from_u64(2);
    assert_eq!(a.pow(&[10]).as_limbs()[0], 17);
}

#[test]
fn legendre_of_qr() {
    // 4 = 2² is a quadratic residue mod 19 → Legendre = 1
    let a = F19::from_u64(4);
    assert_eq!(a.legendre(), 1);
}

#[test]
fn legendre_of_zero() {
    assert_eq!(F19::zero().legendre(), 0);
}

#[test]
fn sqrt_of_qr() {
    // √4 = 2  (mod 19)
    let four = F19::from_u64(4);
    let root = four.sqrt().expect("4 is a QR mod 19");
    assert_eq!((root * root).as_limbs()[0], 4);
}