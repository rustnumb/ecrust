//! Fixed-width secret scalar container.
//!
//! Protocol code should avoid passing variable-length scalar slices for secret
//! material, because the length itself can become observable. This wrapper
//! keeps scalar width fixed at the type level.

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// Fixed-width scalar represented as little-endian `u64` limbs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SecretScalar<const LIMBS: usize> {
    limbs: [u64; LIMBS],
}

impl<const LIMBS: usize> SecretScalar<LIMBS> {
    /// Construct a scalar from its little-endian limb representation.
    pub const fn new(limbs: [u64; LIMBS]) -> Self {
        Self { limbs }
    }

    /// Borrow the underlying little-endian limbs.
    pub fn as_limbs(&self) -> &[u64; LIMBS] {
        &self.limbs
    }

    /// Consume the scalar and return the underlying limbs.
    pub fn to_limbs(self) -> [u64; LIMBS] {
        self.limbs
    }

    /// Total bit width of the scalar container.
    pub const fn bit_len() -> usize {
        64 * LIMBS
    }

    /// Return the bit at absolute big-endian index `i`.
    ///
    /// `i = 0` is the most-significant bit of the highest limb and
    /// `i = bit_len() - 1` is the least-significant bit of limb 0.
    pub fn bit_be(&self, i: usize) -> Choice {
        debug_assert!(i < Self::bit_len());

        let limb_from_msb = i / 64;
        let bit_in_limb = 63 - (i % 64);
        let limb_index = LIMBS - 1 - limb_from_msb;
        Choice::from(((self.limbs[limb_index] >> bit_in_limb) & 1) as u8)
    }
}

impl<const LIMBS: usize> Default for SecretScalar<LIMBS> {
    fn default() -> Self {
        Self { limbs: [0u64; LIMBS] }
    }
}

impl<const LIMBS: usize> ConditionallySelectable for SecretScalar<LIMBS> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let mut limbs = [0u64; LIMBS];
        for i in 0..LIMBS {
            limbs[i] = u64::conditional_select(&a.limbs[i], &b.limbs[i], choice);
        }
        Self { limbs }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        for i in 0..LIMBS {
            self.limbs[i].conditional_assign(&other.limbs[i], choice);
        }
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        for i in 0..LIMBS {
            u64::conditional_swap(&mut a.limbs[i], &mut b.limbs[i], choice);
        }
    }
}

impl<const LIMBS: usize> ConstantTimeEq for SecretScalar<LIMBS> {
    fn ct_eq(&self, other: &Self) -> Choice {
        let mut acc = Choice::from(1u8);
        for i in 0..LIMBS {
            acc = acc & self.limbs[i].ct_eq(&other.limbs[i]);
        }
        acc
    }
}
