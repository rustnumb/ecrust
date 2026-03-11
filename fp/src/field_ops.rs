//! Core trait that every field element in the tower must implement.
//!
//! [`FieldOps`] provides the algebraic contract for additive and
//! multiplicative group operations, together with the utilities
//! (inversion, exponentiation, square root, Frobenius) that
//! algorithms like Vélu's formulas and pairing computations rely on.

use std::ops::{Add, Sub, Mul, Neg};

/// Algebraic operations required of every element in an extension field.
///
/// Implementors live at every level of the tower:
/// [`FpElement`](crate::fp_element::FpElement), `Fp2Element`, `Fp4Element`,
/// `Fp6Element`, and `Fp12Element`.
pub trait FieldOps:
    Sized
    + Clone
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
{
    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    /// The additive identity (zero).
    fn zero() -> Self;

    /// The multiplicative identity (one).
    fn one() -> Self;

    // -----------------------------------------------------------------------
    // Predicates
    // -----------------------------------------------------------------------

    /// Returns `true` iff `self` is the additive identity.
    fn is_zero(&self) -> bool;

    /// Returns `true` iff `self` is the multiplicative identity.
    fn is_one(&self) -> bool;

    // -----------------------------------------------------------------------
    // Arithmetic
    // -----------------------------------------------------------------------

    /// Additive inverse: `−self`.
    fn negate(&self) -> Self;

    /// `self + rhs` (alias kept for use in generic code without the trait bound).
    fn add(&self, rhs: &Self) -> Self;

    /// `self − rhs`.
    fn sub(&self, rhs: &Self) -> Self;

    /// `self × rhs`.
    fn mul(&self, rhs: &Self) -> Self;

    /// `self²` — may be implemented more efficiently than `self.mul(self)`.
    fn square(&self) -> Self;

    /// `self × 2` — cheaper than a full multiplication.
    fn double(&self) -> Self;

    // -----------------------------------------------------------------------
    // Multiplicative inverse / division
    // -----------------------------------------------------------------------

    /// Multiplicative inverse.
    ///
    /// Returns `None` when `self.is_zero()`.
    fn invert(&self) -> Option<Self>;

    /// `self / rhs`.
    ///
    /// Returns `None` when `rhs.is_zero()`.
    fn div(&self, rhs: &Self) -> Option<Self> {
        rhs.invert().map(|inv| self.mul(&inv))
    }

    // -----------------------------------------------------------------------
    // Exponentiation
    // -----------------------------------------------------------------------

    /// `self^exp` using square-and-multiply (little-endian bit order).
    ///
    /// `exp` is supplied as a slice of `u64` limbs, least-significant first.
    fn pow(&self, exp: &[u64]) -> Self {
        let mut result = Self::one();
        let mut base   = self.clone();

        for &limb in exp {
            let mut limb = limb;
            for _ in 0..64 {
                if limb & 1 == 1 {
                    result = result.mul(&base);
                }
                base  = base.square();
                limb >>= 1;
            }
        }
        result
    }

    // -----------------------------------------------------------------------
    // Frobenius endomorphism  φ_p : x ↦ x^p
    // -----------------------------------------------------------------------

    /// Apply the *p*-power Frobenius endomorphism once.
    ///
    /// For `FpElement` this is the identity.
    /// For extension fields it acts as complex conjugation / the Galois action.
    fn frobenius(&self) -> Self;

    /// Apply the Frobenius endomorphism `k` times: `x ↦ x^(p^k)`.
    fn frobenius_pow(&self, k: u32) -> Self {
        let mut result = self.clone();
        for _ in 0..k {
            result = result.frobenius();
        }
        result
    }

    // -----------------------------------------------------------------------
    // Norm and trace over Fp
    // -----------------------------------------------------------------------

    /// Field norm `N_{K/Fp}(self)`.
    ///
    /// For degree-*d* extension: product of all conjugates `x^(p^i)`, `i=0..d-1`.
    fn norm(&self) -> Self;

    /// Field trace `Tr_{K/Fp}(self)`.
    ///
    /// For degree-*d* extension: sum of all conjugates `x^(p^i)`, `i=0..d-1`.
    fn trace(&self) -> Self;

    // -----------------------------------------------------------------------
    // Square root
    // -----------------------------------------------------------------------

    /// Tonelli–Shanks / Cipolla square root in the field.
    ///
    /// Returns `None` if `self` is not a quadratic residue.
    fn sqrt(&self) -> Option<Self>;

    /// Legendre symbol: `1` if QR, `−1` if QNR, `0` if zero.
    fn legendre(&self) -> i8;

    // -----------------------------------------------------------------------
    // Utilities
    // -----------------------------------------------------------------------

    /// Return the characteristic *p* as a little-endian `u64` limb vector.
    fn characteristic() -> Vec<u64>;

    /// Degree of this field over the prime subfield Fp.
    fn degree() -> u32;
}
