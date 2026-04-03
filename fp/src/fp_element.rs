//! Base prime-field element Fp = Z / pZ backed by `crypto-bigint`.

use core::ops::{Add, Mul, Neg, Sub};

use crypto_bigint::{
    modular::{ConstMontyForm, ConstPrimeMontyParams},
    Uint,
};
use subtle::{CtOption, Choice, ConditionallySelectable, ConstantTimeEq};
use crate::field_ops::FieldOps;

/// An element of the prime field Fp = Z/pZ, stored in Montgomery form.
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
    pub fn from_uint(x: Uint<LIMBS>) -> Self {
        Self { value: ConstMontyForm::<MOD, LIMBS>::new(&x) }
    }

    pub fn from_words(words: [u64; LIMBS]) -> Self {
        Self::from_uint(Uint::<LIMBS>::from_words(words))
    }

    pub fn from_limbs(limbs: &[u64]) -> Self {
        assert_eq!(limbs.len(), LIMBS, "wrong number of limbs");
        let mut words = [0u64; LIMBS];
        words.copy_from_slice(limbs);
        Self::from_words(words)
    }

    pub fn from_u64(val: u64) -> Self {
        Self::from_uint(Uint::<LIMBS>::from_u64(val))
    }

    pub fn as_uint(&self) -> Uint<LIMBS> {
        self.value.retrieve()
    }

    pub fn as_limbs(&self) -> [u64; LIMBS] {
        self.value.retrieve().to_words()
    }

    pub fn to_montgomery(&self) -> Uint<LIMBS> {
        self.value.to_montgomery()
    }

    pub fn from_montgomery(mont: Uint<LIMBS>) -> Self {
        Self { value: ConstMontyForm::<MOD, LIMBS>::from_montgomery(mont) }
    }
}



// ---------------------------------------------------------------------------
// CtOption functionalities
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> Default for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn default() -> Self { Self::zero() }
}

impl<MOD, const LIMBS: usize> ConditionallySelectable for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self{ value: ConstMontyForm::conditional_select(&a.value, &b.value, choice) }
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        ConstMontyForm::conditional_swap(&mut a.value, &mut b.value, choice)
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.value.conditional_assign(&other.value, choice)
    }
}

impl<MOD, const LIMBS: usize> ConstantTimeEq for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>
{
    fn ct_eq(&self, other: &Self) -> Choice {
        Choice::from(ConstMontyForm::ct_eq(&self.value, &other.value))
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        Choice::from(ConstMontyForm::ct_ne(&self.value, &other.value))
    }
}


// ---------------------------------------------------------------------------
// Operator overloads
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> Add for FpElement<MOD, LIMBS>
where MOD: ConstPrimeMontyParams<LIMBS>
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self { Self { value: self.value + rhs.value } }
}

impl<MOD, const LIMBS: usize> Sub for FpElement<MOD, LIMBS>
where MOD: ConstPrimeMontyParams<LIMBS>
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self { Self { value: self.value - rhs.value } }
}

impl<MOD, const LIMBS: usize> Mul for FpElement<MOD, LIMBS>
where MOD: ConstPrimeMontyParams<LIMBS>
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self { Self { value: self.value * rhs.value } }
}

impl<MOD, const LIMBS: usize> Neg for FpElement<MOD, LIMBS>
where MOD: ConstPrimeMontyParams<LIMBS>
{
    type Output = Self;
    fn neg(self) -> Self { Self { value: -self.value } }
}

// ---------------------------------------------------------------------------
// FieldOps
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> FieldOps for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn zero() -> Self { Self { value: ConstMontyForm::<MOD, LIMBS>::ZERO } }
    fn one()  -> Self { Self { value: ConstMontyForm::<MOD, LIMBS>::ONE  } }

    fn is_zero(&self) -> Choice { Self::ct_eq(self, &Self::zero()) }
    fn is_one (&self) -> Choice { Self::ct_eq(self, &Self::one()) }

    fn negate(&self) -> Self { Self { value: -self.value } }
    fn add(&self, rhs: &Self) -> Self { Self { value: self.value + rhs.value } }
    fn sub(&self, rhs: &Self) -> Self { Self { value: self.value - rhs.value } }
    fn mul(&self, rhs: &Self) -> Self { Self { value: self.value * rhs.value } }
    fn square(&self) -> Self { Self { value: self.value.square() } }
    fn double(&self) -> Self { Self { value: self.value.double() } }

    fn invert(&self) -> CtOption<Self> {
        let ct_opt_inv = self.value.invert();
        let inv = ct_opt_inv.unwrap();
        let is_invertible = ct_opt_inv.is_some();
        CtOption::new(Self{ value: inv}, is_invertible)
    }

    fn frobenius(&self) -> Self { *self }
    fn norm(&self)      -> Self { *self }
    fn trace(&self)     -> Self { *self }

    fn sqrt(&self) -> Option<Self> {
        self.value.sqrt().into_option().map(|x| Self { value: x })
    }

    fn legendre(&self) -> i8 { i8::from(self.value.jacobi_symbol()) }

    fn characteristic() -> Vec<u64> {
        // We avoid accessing MOD::MODULUS directly because its location in
        // the crypto-bigint trait hierarchy varies across versions.
        //
        // Arithmetic trick: in any field, −1 ≡ p−1 (mod p).
        // Retrieving the Montgomery form of −1 gives p−1 as a plain Uint,
        // so p = (p−1) + 1.  This uses only stable ConstMontyForm constants.
        let minus_one = ConstMontyForm::<MOD, LIMBS>::ZERO - ConstMontyForm::<MOD, LIMBS>::ONE;
        let p_minus_1: Uint<LIMBS> = minus_one.retrieve();
        let p = p_minus_1.wrapping_add(&Uint::<LIMBS>::ONE);
        p.to_words().to_vec()
    }

    fn degree() -> u32 { 1 }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::{const_prime_monty_params, Uint};

    const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2);
    type F19 = FpElement<Fp19Modulus, 1>;

    #[test] fn zero_is_zero()            { assert!(F19::zero().is_zero()); }
    #[test] fn one_is_one()              { assert!(F19::one().is_one().unwrap()); }
    #[test] fn degree_is_one()           { assert_eq!(F19::degree(), 1); }
    #[test] fn characteristic_is_19()   { assert_eq!(F19::characteristic(), vec![19u64]); }

    #[test] fn add_mod_p() {
        // 17 + 5 = 22 ≡ 3 (mod 19)
        assert_eq!((F19::from_u64(17) + F19::from_u64(5)).as_limbs()[0], 3);
    }

    #[test] fn sub_mod_p() {
        // 3 − 7 = −4 ≡ 15 (mod 19)
        assert_eq!((F19::from_u64(3) - F19::from_u64(7)).as_limbs()[0], 15);
    }

    #[test] fn mul_mod_p() {
        // 7 × 8 = 56 ≡ 18 (mod 19)
        assert_eq!((F19::from_u64(7) * F19::from_u64(8)).as_limbs()[0], 18);
    }

    #[test] fn neg_mod_p() {
        // −3 ≡ 16 (mod 19)
        assert_eq!((-F19::from_u64(3)).as_limbs()[0], 16);
    }

    #[test] fn inv_works() {
        let a = F19::from_u64(7);
        assert_eq!((a * a.invert().unwrap()).as_limbs()[0], 1);
    }

    #[test] fn inv_zero_is_none() { assert!(F19::zero().invert().is_none()); }

    #[test] fn pow_works() {
        // 2^10 = 1024 ≡ 17 (mod 19)
        use crate::field_ops::FieldOps;
        assert_eq!(F19::from_u64(2).pow(&[10]).as_limbs()[0], 17);
    }

    #[test] fn sqrt_of_qr() {
        // √4 squared must give back 4
        let four = F19::from_u64(4);
        let root = four.sqrt().expect("4 is a QR mod 19");
        assert_eq!((root * root).as_limbs()[0], 4);
    }
}