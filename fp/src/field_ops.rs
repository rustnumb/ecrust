//! Core trait that every field element in the tower must implement.

use std::ops::{Add, Mul, Neg, Sub};
use subtle::{Choice, ConditionallySelectable, CtOption};

/// Trait for generating cryptographically secure random field elements.
///
/// Separated from [`FieldOps`] so that downstream code that doesn't
/// need randomness is free of the `rand` dependency.
pub trait FieldRandom: Sized {
    /// Sample a uniformly random element using a cryptographic RNG.
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self;
}

pub trait FieldOps:
    Sized
    + Clone
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Default
    + ConditionallySelectable
// const time impl will have the following trait too
// + ConditionallySelectable
{
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> Choice;
    fn is_one(&self) -> Choice;
    fn negate(&self) -> Self;
    fn add(&self, rhs: &Self) -> Self;
    fn sub(&self, rhs: &Self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;
    fn square(&self) -> Self;
    fn double(&self) -> Self;
    fn invert(&self) -> CtOption<Self>;

    fn div(&self, rhs: &Self) -> CtOption<Self> {
        rhs.invert().map(|inv| self.mul(&inv))
    }

    /// `self^exp` using square-and multiply (litte-endian bit order)
    ///
    /// It is constant time for fixed `exp`
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
    /// WARNING: Only constant time if the number of limbs of exp is
    /// constant
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

    /// Returns a squareroot if it exists
    ///
    /// Returns a squareroof of `self` if it exists in the finite
    /// field FpM. The return type is Ctoption<Self> and it is not
    /// none if and only if the squareroot exists. This is an
    /// implementation fo the Tonelli--Shanks algorithm in the
    /// multiplicative group FpM*
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of FpM (type: self)
    ///
    /// # Returns
    ///
    /// The square root of `self` (type: CtOption<Self>)
    fn sqrt(&self) -> CtOption<Self>;

    /// Computes the inverse and square root of `self`
    ///
    /// Computes simulaineously the inverse and square root of `self`.
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of FpM (type: self)
    ///
    /// # Returns
    ///
    /// The inverse and square root of `self`. The former is none if
    /// and only if nonzero and the latter is not none if and only if
    /// there exists a squareroot in FpM
    /// (type: (CtOption<Self>, CtOption<self>))
    fn inverse_and_sqrt(&self) -> (CtOption<Self>, CtOption<Self>) {
        (self.invert(), self.sqrt())
    }

    /// Computes the square root the inverse of `self`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of FpM (type: self)
    ///
    /// # Returns
    ///
    /// The square root of the inverse of `self`. The former is not
    /// none if and only if it is both nonzero there exists a
    /// squareroot in FpM (type: CtOption<self>)
    fn inv_sqrt(&self) -> CtOption<Self> {
        self.sqrt().and_then(|s| s.invert())
    }

    /// Computes the inverse of `self` and square root of `rhs`
    ///
    /// Computes simulaineously the inverse of `self` and square root
    /// of `rhs`.
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of FpM (type: self)
    /// * `rhs` - Element of FpM (type: &Self)
    ///
    /// # Returns
    ///
    /// The inverse of `self` and square root fo `rhs`. Theq former is
    /// none if and only if `self` is nonzero and the latter is not
    /// none if and only if there exists a squareroot of `rhs` in FpM
    /// (type: (CtOption<Self>, CtOption<self>))
    fn invertme_sqrtother(&self, rhs: &Self) -> (CtOption<Self>, CtOption<Self>) {
        (self.invert(), rhs.sqrt())
    }

    /// Computes the squareroot of a ratio `self/rhs`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of FpM (type: self)
    /// * `rhs` - Element of FpM (type: &Self)
    ///
    /// # Returns
    ///
    /// The squareroot of the ratio `self/rhs` is not none if and only
    /// if `rhs` is invertible and the ratio has an FpM squareroot
    /// (type: (CtOption<Self>, CtOption<self>))
    fn sqrt_ratio(&self, rhs: &Self) -> CtOption<Self> {
        rhs.invert().and_then(|inv_rhs| self.mul(&inv_rhs).sqrt())
    }


    fn legendre(&self) -> i8;
    fn characteristic() -> Vec<u64>;
    fn degree() -> u32;
}
