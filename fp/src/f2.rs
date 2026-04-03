//! Base binary field element F2 = Z / 2Z

use core::ops::{Add, Mul, Neg, Sub};
use crypto_bigint::Uint;
use subtle::{Choice, CtOption, ConditionallySelectable, ConstantTimeEq};
use crate::field_ops::FieldOps;


#[derive(Debug, Eq, PartialEq, Copy, Clone)]
pub struct F2Element {
    pub(crate) value: Uint<1>
}



// ---------------------------------------------------------------------------
// Constructors / accessors
// ---------------------------------------------------------------------------

impl F2Element {
    pub const ZERO: Self = Self { value: Uint::<1>::ZERO };
    pub const ONE:  Self = Self { value: Uint::<1>::ONE };

    fn new(x: Uint<1>) -> Self { Self { value: x & Uint::<1>::ONE } }

    pub fn from_u64(x : u64) -> Self{
        Self::new(Uint::from(x & 1))
    }

    pub fn value(&self) -> Uint<1> {
        self.value
    }

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

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        Uint::<1>::conditional_swap(&mut a.value, &mut b.value, choice)
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.value.conditional_assign(&other.value, choice)
    }
}

impl ConstantTimeEq for F2Element {
    fn ct_eq(&self, other: &Self) -> Choice {
        Choice::from(Uint::<1>::ct_eq(&self.value, &other.value))
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        Choice::from(Uint::<1>::ct_ne(&self.value, &other.value))
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
    fn zero() -> Self { Self::ZERO }
    fn one() -> Self { Self::ONE }

    fn is_zero(&self) -> Choice { Self::ct_eq(self, &Self::zero()) }
    fn is_one(&self) -> Choice { Self::ct_eq(self, &Self::one()) }
    fn negate(&self) -> Self { *self }
    fn add(&self, rhs: &Self) -> Self { Self::new(self.value ^ rhs.value) }
    fn sub(&self, rhs: &Self) -> Self { Self::new(self.value ^ rhs.value) }
    fn mul(&self, rhs: &Self) -> Self { Self::new(self.value & rhs.value) }
    fn square(&self) -> Self { *self }      // x = x^2 for every x in F_2
    fn double(&self) -> Self { Self::ZERO }
    fn invert(&self) -> CtOption<Self> {
        let is_invertible = !self.is_zero();
        CtOption::new(*self, is_invertible)
    }

    fn frobenius(&self) -> Self { *self }
    fn norm(&self)      -> Self { *self }
    fn trace(&self)     -> Self { *self }

    fn sqrt(&self) -> Option<Self> { Some(*self)}

    fn legendre(&self) -> i8 {
        self.as_u8() as i8
    }

    fn characteristic() -> Vec<u64> { vec![2] }
    fn degree() -> u32 { 1 }
}
