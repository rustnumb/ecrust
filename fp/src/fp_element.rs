//! Base prime-field element $\mathbb{F}_p = \mathbb{Z} / p\mathbb{Z}$
//! backed by `crypto-bigint`.
//!
//! # Overview
//!
//! Given an odd prime $p$ supplies $\mathbb{F}_p$ the finite field with
//! $p$ elements.
//!
//! This module provides a single generic type [`FpElement<MOD, const
//! LIMBS: usize>`](self::fp_element::FpElement).
//!
//! # Example
//!
//! ```
//! use crypto_bigint::{const_prime_monty_params, Uint};
//! use fp::fp_element::FpElement;
//! use fp::field_ops::FieldOps;
//!
//! /*
//! We will set up the field F_19 so note that 0000000000000013 is
//! equal to 19 in Hexadecimal and that 2 is a generator (primitive
//! root) of (Z/pZ)^*
//! */
//!
//! const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2);
//! type F19 = FpElement<Fp19Modulus, 1>;
//!
//! /* Some standard tests */
//! let a = F19::from_u64(17);
//! let b = F19::from_u64(5);
//! assert_eq!((&a + &b).as_limbs()[0], 3);
//! assert_eq!((&a - &b).as_limbs()[0], 12);
//! assert_eq!((&a * &b).as_limbs()[0], 9);
//! assert_eq!((a.invert().unwrap()).as_limbs()[0], 9);
//! ```

use core::ops::{Add, Mul, Neg, Sub};
use std::fmt;

use crate::field_ops::{FieldFromRepr, FieldOps, FieldRandom};
use crypto_bigint::{
    modular::{ConstMontyForm, ConstPrimeMontyParams},
    NonZero, RandomMod, Uint,
};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// An element of the prime field $\mathbb{F}_p =
/// \mathbb{Z}/p\mathbb{Z}$, stored in Montgomery form.
///
/// # Note
///
/// The internal value uses `crypto-bigint`'s [`ConstMontyForm`], so arithmetic
/// is performed in Montgomery representation while the public constructors and
/// accessors accept and return canonical integers.
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
    /// Goes from an `Uint<LIMBS>` to an element of $\mathbb{F}_p$
    ///
    /// # Arguments
    ///
    /// * `x` - An integer (type: `Uint<LIMBS>`)
    ///
    /// # Returns
    ///
    /// Element of $\mathbb{F}_p$ (type: `Self`)
    pub fn from_uint(x: Uint<LIMBS>) -> Self {
        Self {
            value: ConstMontyForm::<MOD, LIMBS>::new(&x),
        }
    }

    /// Goes from an `[u64; LIMBS]` to an element of $\mathbb{F}_p$
    ///
    /// # Arguments
    ///
    /// * `words` - A vec of `u64` (type: `[u64; LIMBS]`)
    ///
    /// # Returns
    ///
    /// Element of $\mathbb{F}_p$ (type: `Self`)
    pub fn from_words(words: [u64; LIMBS]) -> Self {
        Self::from_uint(Uint::<LIMBS>::from_words(words))
    }

    /// Goes from an `&[u64]` to an element of $\mathbb{F}_p$
    ///
    /// # Arguments
    ///
    /// * `limbs` - A vec of `u64` (type: `&[u64]`)
    ///
    /// # Returns
    ///
    /// Element of $\mathbb{F}_p$ (type: `Self`)
    pub fn from_limbs(limbs: &[u64]) -> Self {
        assert_eq!(limbs.len(), LIMBS, "wrong number of limbs");
        let mut words = [0u64; LIMBS];
        words.copy_from_slice(limbs);
        Self::from_words(words)
    }

    /// Goes from an `u64` to an element of $\mathbb{F}_p$
    ///
    /// # Arguments
    ///
    /// * `val` - A `u64` int (type: `u64`)
    ///
    /// # Returns
    ///
    /// Element of $\mathbb{F}_p$ (type: `Self`)
    pub fn from_u64(val: u64) -> Self {
        Self::from_uint(Uint::<LIMBS>::from_u64(val))
    }

    /// Goes from an element of $\mathbb{F}_p$ to the corresponding `Uint`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_p$ (type: `Self`)
    ///
    /// # Returns
    ///
    /// The corresponding `Uint` (type: `Uint<LIMBS>`)
    pub fn as_uint(&self) -> Uint<LIMBS> {
        self.value.retrieve()
    }

    /// Goes from an element of $\mathbb{F}_p$ to the limbs
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_p$ (type: `Self`)
    ///
    /// # Returns
    ///
    /// The limbs giving the integer (type: `[u64; LIMBS]`)
    pub fn as_limbs(&self) -> [u64; LIMBS] {
        self.value.retrieve().to_words()
    }

    /// Gives a montgomery `Uint<LIMBS>` from an $\mathbb{F}_p$
    /// element
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_p$ (type: `Self`)
    ///
    /// # Returns
    ///
    /// Integer in montgomery from (type: `Uint<LIMBS>`)
    pub fn to_montgomery(&self) -> Uint<LIMBS> {
        self.value.to_montgomery()
    }

    /// Gives an element of $\mathbb{F}_p$ from the Montgomery
    /// representation.
    ///
    /// # Arguments
    ///
    /// * `mont` - The input montgomery value (type: `Uint<LIMBS>`)
    ///
    /// # Returns
    ///
    /// Corresponding element of $\mathbb{F}_p$ (type: `Self`)
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

// Manual `Hash` via the canonical (non-Montgomery) representative, so that
// `a == b` ⇒ `hash(a) == hash(b)` independently of any internal Montgomery
// encoding.  Matches the semantics of the derived `PartialEq` above, which
// ultimately compares `ConstMontyForm` values (themselves bijective with
// their canonical form).
impl<MOD, const LIMBS: usize> core::hash::Hash for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.as_uint().hash(state);
    }
}

// ---------------------------------------------------------------------------
// Operator overloads
// ---------------------------------------------------------------------------

impl<'a, 'b, MOD, const LIMBS: usize> Add<&'b FpElement<MOD, LIMBS>> for &'a FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = FpElement<MOD, LIMBS>;

    fn add(self, rhs: &'b FpElement<MOD, LIMBS>) -> Self::Output {
        <FpElement<MOD, LIMBS> as FieldOps>::add(self, rhs)
    }
}

impl<'a, 'b, MOD, const LIMBS: usize> Sub<&'b FpElement<MOD, LIMBS>> for &'a FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = FpElement<MOD, LIMBS>;

    fn sub(self, rhs: &'b FpElement<MOD, LIMBS>) -> Self::Output {
        <FpElement<MOD, LIMBS> as FieldOps>::sub(self, rhs)
    }
}

impl<'a, 'b, MOD, const LIMBS: usize> Mul<&'b FpElement<MOD, LIMBS>> for &'a FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = FpElement<MOD, LIMBS>;

    fn mul(self, rhs: &'b FpElement<MOD, LIMBS>) -> Self::Output {
        <FpElement<MOD, LIMBS> as FieldOps>::mul(self, rhs)
    }
}

impl<'a, MOD, const LIMBS: usize> Neg for &'a FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Output = FpElement<MOD, LIMBS>;

    fn neg(self) -> Self::Output {
        <FpElement<MOD, LIMBS> as FieldOps>::negate(self)
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

    fn from_u64(x: u64) -> Self {
        Self::from_u64(x)
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

    fn legendre(&self) -> i8 {
        i8::from(self.value.jacobi_symbol())
    }

    fn characteristic() -> Vec<u64> {
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
        if LIMBS == 1 {
            write!(f, "{}", val.to_words()[0])
        } else {
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
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self {
        let minus_one = ConstMontyForm::<MOD, LIMBS>::ZERO - ConstMontyForm::<MOD, LIMBS>::ONE;
        let p_minus_1: Uint<LIMBS> = minus_one.retrieve();
        let p = p_minus_1.wrapping_add(&Uint::<LIMBS>::ONE);
        let modulus = NonZero::new(p).expect("prime modulus must be nonzero");
        let val = Uint::<LIMBS>::random_mod_vartime(rng, &modulus);
        Self::from_uint(val)
    }
}

impl<MOD, const LIMBS: usize> FieldFromRepr for FpElement<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    type Repr = Uint<LIMBS>;

    fn from_repr(x: Self::Repr) -> Self {
        Self::from_uint(x)
    }
}
