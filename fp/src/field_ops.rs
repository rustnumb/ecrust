//! Core trait that every field element in the tower must implement.

use std::ops::{Add, Mul, Neg, Sub};
// use subtle::ConditionallySelectable;

pub trait FieldOps:
    Sized
    + Clone
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
// const time impl will have the following trait too
// + ConditionallySelectable
{
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn negate(&self) -> Self;
    fn add(&self, rhs: &Self) -> Self;
    fn sub(&self, rhs: &Self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;
    fn square(&self) -> Self;
    fn double(&self) -> Self;
    fn invert(&self) -> Option<Self>;

    fn div(&self, rhs: &Self) -> Option<Self> {
        rhs.invert().map(|inv| self.mul(&inv))
    }

    /// `self^exp` using square-and multiply (litte-endian bit order)
    ///
    /// # Arguments
    ///
    /// * `&self` - Finite field element (type: self)
    /// * `exp` - Exponent (type: &[u64])
    ///
    /// # Returns
    ///
    /// `&self^exp` (type: Self)
    ///
    /// # Todo
    ///
    /// Implement constant time version.
    ///
    /// # Why `<Self as FieldOps>::mul` instead of `result.mul(&base)`
    ///
    /// `FieldOps` requires `Mul<Output = Self>` as a supertrait, so `Self`
    /// exposes **two** methods named `mul`:
    ///
    ///   - `<Self as Mul>::mul(self, rhs: Self) -> Self`   ← operator, takes by value
    ///   - `<Self as FieldOps>::mul(&self, rhs: &Self) -> Self` ← ours, takes by ref
    ///
    /// Writing `result.mul(&base)` triggers method resolution, which picks
    /// `Mul::mul` (the operator) because it was declared first in the supertrait
    /// list. `Mul::mul` expects `Self`, not `&Self` → E0308.
    ///
    /// Fully-qualified syntax `<Self as FieldOps>::mul(...)` bypasses method
    /// resolution entirely and calls exactly the trait method we want.
    fn pow_vartime(&self, exp: &[u64]) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();

        for &limb in exp {
            let mut limb = limb;
            for _ in 0..64 {
                if limb & 1 == 1 {
                    result = <Self as FieldOps>::mul(&result, &base);
                }
                base = <Self as FieldOps>::square(&base);
                limb >>= 1;
            }
        }
        result
    }

    /// `self^pow` in constant time using a Montgomery ladder
    ///
    /// Uses a Montgomery ladder to compute `self^exp`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_p$ (type: self)
    /// * `exp` - Exponent (type: &[u64])
    ///
    /// # Returns
    ///
    /// The value `self^pow` (type: Self)
    ///
    /// # Todo
    ///
    /// Use `subtle` and `conditional_swap` to make true constant time
    fn pow(&self, exp: &[u64]) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();

        for &limb in exp.iter().rev() {
            let mut limb = limb.reverse_bits();
            for _ in 0..64 {
                let mybit = limb & 1;
                // TODO: true constant time implementation will have the below
                // conditional_swap(&base, &result, ((limb & 1) as u8).into());
                if mybit == 1 {
                    (result, base) = (base, result);
                }
                base = <Self as FieldOps>::mul(&result, &base);
                result = <Self as FieldOps>::square(&result);
                // conditional_swap(&base, &result, ((limb & 1) as u8).into());
                if mybit == 1 {
                    (result, base) = (base, result);
                }
                limb >>= 1;
            }
        }
        result
    }

    fn frobenius(&self) -> Self;

    fn frobenius_pow(&self, k: u32) -> Self {
        let mut result = self.clone();
        for _ in 0..k {
            result = result.frobenius();
        }
        result
    }

    fn norm(&self) -> Self;
    fn trace(&self) -> Self;
    fn sqrt(&self) -> Option<Self>;
    fn legendre(&self) -> i8;
    fn characteristic() -> Vec<u64>;
    fn degree() -> u32;
}
