//! Base prime-field element Fp = Z / pZ backed by `crypto-bigint`.

use core::ops::{Add, Mul, Neg, Sub};

use crypto_bigint::{
    modular::{ConstMontyForm, ConstPrimeMontyParams},
    Uint,
};

use crate::field_ops::FieldOps;

/// An element of the prime field Fp = Z/pZ.
///
/// Internally stored in Montgomery form via `ConstMontyForm`.
/// This version assumes a compile-time modulus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FpElement<MOD, const LIMBS: usize>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    pub(crate) value: ConstMontyForm<MOD, LIMBS>,
}

// ---------------------------------------------------------------------------
// Constructors / accessors
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    /// Construct from a reduced or unreduced integer; reduction mod p is applied.
    pub fn from_uint(x: Uint<LIMBS>) -> Self {
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::new(&x),
        }
    }

    /// Construct from little-endian machine words.
    ///
    /// On 64-bit targets, these are `u64` words.
    pub fn from_words(words: [u64; LIMBS]) -> Self {
        Self::from_uint(Uint::<LIMBS>::from_words(words))
    }

    /// Construct from a slice of little-endian words.
    pub fn from_limbs(limbs: &[u64]) -> Self {
        assert_eq!(limbs.len(), LIMBS, "wrong number of limbs");
        let mut words = [0u64; LIMBS];
        words.copy_from_slice(limbs);
        Self::from_words(words)
    }

    /// Construct from a single `u64`.
    pub fn from_u64(val: u64) -> Self {
        Self::from_uint(Uint::<LIMBS>::from_u64(val))
    }

    /// Return the canonical representative in `[0, p)`.
    pub fn as_uint(&self) -> Uint<LIMBS> {
        self.value.retrieve()
    }

    /// Return canonical little-endian words.
    pub fn as_limbs(&self) -> [u64; LIMBS] {
        self.value.retrieve().to_words()
    }

    /// Return the raw Montgomery representation.
    pub fn to_montgomery(&self) -> Uint<LIMBS> {
        self.value.to_montgomery()
    }

    /// Build directly from a Montgomery-form integer.
    pub fn from_montgomery(mont: Uint<LIMBS>) -> Self {
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::from_montgomery(mont),
        }
    }
}

// ---------------------------------------------------------------------------
// Operator overloads
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> Add for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
        }
    }
}

impl<MOD, const LIMBS: usize> Sub for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            value: self.value - rhs.value,
        }
    }
}

impl<MOD, const LIMBS: usize> Mul for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
        }
    }
}

impl<MOD, const LIMBS: usize> Neg for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = Self;

    fn neg(self) -> Self {
        Self { value: -self.value }
    }
}

// ---------------------------------------------------------------------------
// FieldOps implementation
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> FieldOps for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn zero() -> Self {
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::ZERO,
        }
    }

    fn one() -> Self {
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::ONE,
        }
    }

    fn is_zero(&self) -> bool {
        self.value == ConstMontyForm::<MOD, LIMBS>::ZERO
    }

    fn is_one(&self) -> bool {
        self.value == ConstMontyForm::<MOD, LIMBS>::ONE
    }

    fn negate(&self) -> Self {
        Self { value: -self.value }
    }

    fn add(&self, rhs: &Self) -> Self {
        Self {
            value: self.value + rhs.value,
        }
    }

    fn sub(&self, rhs: &Self) -> Self {
        Self {
            value: self.value - rhs.value,
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self {
            value: self.value * rhs.value,
        }
    }

    fn square(&self) -> Self {
        Self {
            value: self.value.square(),
        }
    }

    fn double(&self) -> Self {
        Self {
            value: self.value.double(),
        }
    }

    fn invert(&self) -> Option<Self> {
        self.value
            .invert()
            .into_option()
            .map(|x| Self { value: x })
    }

    fn frobenius(&self) -> Self {
        *self
    }

    fn norm(&self) -> Self {
        *self
    }

    fn trace(&self) -> Self {
        *self
    }

    fn sqrt(&self) -> Option<Self> {
        self.value
            .sqrt()
            .into_option()
            .map(|x| Self { value: x })
    }

    fn legendre(&self) -> i8 {
        i8::from(self.value.jacobi_symbol())
    }

    fn characteristic() -> Vec<u64> {
        MOD::MODULUS.as_ref().as_words().to_vec()
    }

    fn degree() -> u32 {
        1
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::{const_prime_monty_params, Uint};

    // p = 19
    const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2, "Fp mod 19");
    type F19 = FpElement<Fp19Modulus, 1>;

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
        let a = F19::from_u64(17);
        let b = F19::from_u64(5);
        assert_eq!((a + b).as_limbs()[0], 3);
    }

    #[test]
    fn mul_mod_p() {
        let a = F19::from_u64(7);
        let b = F19::from_u64(8);
        assert_eq!((a * b).as_limbs()[0], 18);
    }

    #[test]
    fn inv_works() {
        let a = F19::from_u64(7);
        let inv = a.invert().unwrap();
        assert_eq!((a * inv).as_limbs()[0], 1);
    }
}