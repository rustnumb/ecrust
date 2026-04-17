//! Base binary field element $\mathbb{F}_2 = \mathbb{Z} /
//! 2\mathbb{Z}$

use crate::field_ops::FieldFromRepr;
use crate::field_ops::{FieldOps, FieldRandom};
use core::ops::{Add, Mul, Neg, Sub};
use crypto_bigint::Uint;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// Element of the finite field $\mathbb{F}_2$
#[derive(Debug, Eq, PartialEq, Copy, Clone)]
pub struct F2Element {
    pub(crate) value: Uint<1>,
}

// ---------------------------------------------------------------------------
// Constructors / accessors
// ---------------------------------------------------------------------------

impl F2Element {
    /// The constant zero
    pub const ZERO: Self = Self {
        value: Uint::<1>::ZERO,
    };

    /// The constant one
    pub const ONE: Self = Self {
        value: Uint::<1>::ONE,
    };

    /// Create a new element of $\mathbb{F}_2$
    fn new(x: Uint<1>) -> Self {
        Self {
            value: x & Uint::<1>::ONE,
        }
    }

    /// Create a new element of $\mathbb{F}_2$ from a `u64`
    pub fn from_u64(x: u64) -> Self {
        Self::new(Uint::from(x & 1))
    }

    /// Get the `Uint<1>` from an element of $\mathbb{F}_2$
    pub fn value(&self) -> Uint<1> {
        self.value
    }

    /// Get the `u8` from an element of $\mathbb{F}_2$
    pub fn as_u8(&self) -> u8 {
        self.value.to_words()[0] as u8
    }
}

// ---------------------------------------------------------------------------
// CtOption functionalities
// ---------------------------------------------------------------------------

impl Default for F2Element {
    fn default() -> Self {
        Self::zero()
    }
}

impl ConditionallySelectable for F2Element {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self::new(Uint::<1>::conditional_select(&a.value, &b.value, choice))
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.value.conditional_assign(&other.value, choice)
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        Uint::<1>::conditional_swap(&mut a.value, &mut b.value, choice)
    }
}

impl ConstantTimeEq for F2Element {
    fn ct_eq(&self, other: &Self) -> Choice {
        Uint::<1>::ct_eq(&self.value, &other.value)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        Uint::<1>::ct_ne(&self.value, &other.value)
    }
}

// ---------------------------------------------------------------------------
// Operator overloads
// ---------------------------------------------------------------------------

impl Add for F2Element {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new(self.value ^ rhs.value)
    }
}

impl Sub for F2Element {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        //subtraction = addition in F_2
        Self::new(self.value ^ rhs.value)
    }
}

impl Mul for F2Element {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::new(self.value & rhs.value)
    }
}

impl Neg for F2Element {
    type Output = Self;

    fn neg(self) -> Self {
        // x = -x for any x in F_2
        self
    }
}

// ---------------------------------------------------------------------------
// FieldOps implementation
// ---------------------------------------------------------------------------

impl FieldOps for F2Element {
    fn zero() -> Self {
        Self::ZERO
    }
    fn one() -> Self {
        Self::ONE
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
        *self
    }
    fn add(&self, rhs: &Self) -> Self {
        Self::new(self.value ^ rhs.value)
    }
    fn sub(&self, rhs: &Self) -> Self {
        Self::new(self.value ^ rhs.value)
    }
    fn mul(&self, rhs: &Self) -> Self {
        Self::new(self.value & rhs.value)
    }
    fn square(&self) -> Self {
        *self
    } // x = x^2 for every x in F_2
    fn double(&self) -> Self {
        Self::ZERO
    }
    fn invert(&self) -> CtOption<Self> {
        let is_invertible = !self.is_zero();
        CtOption::new(*self, is_invertible)
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
        let is_non_zero = !self.is_zero();
        CtOption::new(*self, is_non_zero)
    }

    fn legendre(&self) -> i8 {
        self.as_u8() as i8
    }

    fn characteristic() -> Vec<u64> {
        vec![2]
    }
    fn degree() -> u32 {
        1
    }
}

// ---------------------------------------------------------------------------
// Pretty Display
// ---------------------------------------------------------------------------

impl std::fmt::Display for F2Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_u8())
    }
}

// ---------------------------------------------------------------------------
// Cryptographically secure random sampling
// ---------------------------------------------------------------------------

impl FieldRandom for F2Element {
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self {
        let bit = (rng.next_u32() & 1) as u64;
        Self::from_u64(bit)
    }
}

impl FieldFromRepr for F2Element {
    type Repr = Uint<1>;

    fn from_repr(x: Self::Repr) -> Self {
        Self::new(x)
    }
}
