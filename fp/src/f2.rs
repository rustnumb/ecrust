//! Base binary field element F2 = Z / 2Z

use core::ops::{Add, Mul, Neg, Sub};
use crypto_bigint::Uint;
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

    fn is_zero(&self) -> bool { *self == Self::ZERO }
    fn is_one(&self) -> bool { *self == Self::ONE }
    fn negate(&self) -> Self { *self }
    fn add(&self, rhs: &Self) -> Self { Self::new(self.value ^ rhs.value) }
    fn sub(&self, rhs: &Self) -> Self { Self::new(self.value ^ rhs.value) }
    fn mul(&self, rhs: &Self) -> Self { Self::new(self.value & rhs.value) }
    fn square(&self) -> Self { *self }      // x = x^2 for every x in F_2
    fn double(&self) -> Self { Self::ZERO }
    fn invert(&self) -> Option<Self> { if self.is_zero() { None } else { Some(*self) } }

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


// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test] fn zero_is_zero() { assert!(F2Element::ZERO.is_zero()); }
    #[test] fn one_is_one() { assert!(F2Element::ONE.is_one()); }
    #[test] fn degree_is_one() { assert_eq!(F2Element::degree(), 1); }
    #[test] fn characteristic_is_2() { assert_eq!(F2Element::characteristic(), vec![2u64]); }
    #[test] fn add_mod_2() {
        // 0 + 1 ≡ 1 (mod 2)
        // 1 + 1 ≡ 0 (mod 2)
        assert_eq!(F2Element::from_u64(0) + F2Element::from_u64(1), F2Element::ONE);
        assert_eq!(F2Element::from_u64(1) + F2Element::from_u64(1), F2Element::ZERO);
    }
    #[test] fn mul_mod_2() {
        // 0 * 0 ≡ 0 (mod 2)
        // 0 * 1 ≡ 0 (mod 2)
        // 1 * 1 ≡ 1 (mod 2)
        assert_eq!(F2Element::from_u64(0) * F2Element::from_u64(0), F2Element::ZERO);
        assert_eq!(F2Element::from_u64(0) * F2Element::from_u64(1), F2Element::ZERO);
        assert_eq!(F2Element::from_u64(1) * F2Element::from_u64(1), F2Element::ONE);
    }
    #[test] fn neg_mod_2()   { assert_eq!(F2Element::ONE.negate(), F2Element::ONE);}
    #[test] fn inv_mod_2()   { assert_eq!(F2Element::ONE.invert().unwrap(), F2Element::ONE);}

}