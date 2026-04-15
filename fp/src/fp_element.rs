//! Base prime-field element Fp = Z / pZ backed by `crypto-bigint`.

use core::ops::{Add, Mul, Neg, Sub};
use std::fmt;

use crate::field_ops::{FieldOps, FieldRandom};
use crypto_bigint::{
    modular::{ConstMontyForm, ConstPrimeMontyParams},
    NonZero, RandomMod, Uint,
};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

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
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::new(&x),
        }
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
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::from_montgomery(mont),
        }
    }
}

// ---------------------------------------------------------------------------
// CtOption functionalities
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> Default for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn default() -> Self {
        Self::zero()
    }
}

impl<MOD, const LIMBS: usize> ConditionallySelectable for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            value: ConstMontyForm::conditional_select(&a.value, &b.value, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.value.conditional_assign(&other.value, choice)
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        ConstMontyForm::conditional_swap(&mut a.value, &mut b.value, choice)
    }
}

impl<MOD, const LIMBS: usize> ConstantTimeEq for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        ConstMontyForm::ct_eq(&self.value, &other.value)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        ConstMontyForm::ct_ne(&self.value, &other.value)
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
// FieldOps
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

    fn is_zero(&self) -> Choice {
        Self::ct_eq(self, &Self::zero())
    }
    fn is_one(&self) -> Choice {
        Self::ct_eq(self, &Self::one())
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

    fn invert(&self) -> CtOption<Self> {
        self.value.invert().map(|inv| Self { value: inv }).into()
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

    fn sqrt(&self) -> CtOption<Self> {
        self.value.sqrt().map(|sqrt| Self { value: sqrt }).into()
    }

    // TODO: Implement fast version
    // fn inv_sqrt(&self) -> CtOption<Self> {
    //
    // }

    fn legendre(&self) -> i8 {
        i8::from(self.value.jacobi_symbol())
    }

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

    fn degree() -> u32 {
        1
    }
}

// ---------------------------------------------------------------------------
// Pretty Display
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> fmt::Display for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let val = self.as_uint();
        // For single-limb fields, show a compact decimal representation.
        if LIMBS == 1 {
            write!(f, "{}", val.to_words()[0])
        } else {
            // Multi-limb: show lower-case hex with 0x prefix.
            write!(f, "0x{val:x}")
        }
    }
}

// ---------------------------------------------------------------------------
// Cryptographically secure random sampling
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize> FieldRandom for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    /// Sample a uniformly random element of Fp using a CSPRNG.
    ///
    /// Internally uses `crypto_bigint::RandomMod::random_mod_vartime` to get a value in `[0, p)`.
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self {
        // Recover the modulus p the same way `characteristic()` does.
        let minus_one = ConstMontyForm::<MOD, LIMBS>::ZERO - ConstMontyForm::<MOD, LIMBS>::ONE;
        let p_minus_1: Uint<LIMBS> = minus_one.retrieve();
        let p = p_minus_1.wrapping_add(&Uint::<LIMBS>::ONE);

        // `random_mod` returns a uniform value in [0, p).
        let modulus = NonZero::new(p).expect("prime modulus must be nonzero");
        let val = Uint::<LIMBS>::random_mod_vartime(rng, &modulus);
        Self::from_uint(val)
    }
}
