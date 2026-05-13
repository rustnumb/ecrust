//! Concrete field adapter for the Jacobi vector tests.

use crypto_bigint::{Uint, const_prime_monty_params};
use fp::fp_element::FpElement;

// Use the same pattern as your Weierstrass tests, but over GF(101).
// Replace the last Montgomery parameter with the one used in your project
// for p = 101 if needed.
const_prime_monty_params!(Fp101Mod, Uint<1>, "0000000000000065", 2);

pub type F = FpElement<Fp101Mod, 1>;

pub fn f(n: i64) -> F {
    let p = 101i64;
    F::from_u64(n.rem_euclid(p) as u64)
}
