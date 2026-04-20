//! Core trait that every field element in the tower must implement.

use subtle::{Choice, ConditionallySelectable, CtOption};

/// Trait for generating cryptographically secure random field elements.
///
/// Separated from [`FieldOps`] so that downstream code that doesn't
/// need randomness is free of the `rand` dependency.
pub trait FieldRandom: Sized {
    /// Sample a uniformly random element using a cryptographic RNG.
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self;
}

/// Core arithmetic interface for field elements used throughout the library.
///
/// This trait abstracts the operations needed by the prime-field layer,
/// extension fields, and higher-level elliptic-curve code.
///
/// It combines:
/// - basic ring-style operations (`+`, `-`, `*`, unary `-`);
/// - distinguished constants `zero()` and `one()`;
/// - predicates such as `is_zero()` and `is_one()`;
/// - field-specific operations such as inversion, Frobenius, norm, trace,
///   square roots, and Legendre-symbol-style square testing;
/// - constant-time conditional selection via `subtle::ConditionallySelectable`.
///
/// The trait is intentionally separate from [`FieldRandom`], so downstream code
/// that only needs deterministic arithmetic does not need to depend on `rand`.
///
/// Scalars used by exponentiation methods are encoded as little-endian `u64`
/// limbs, matching the convention used elsewhere in the library.
pub trait FieldOps:
    Sized
    + Clone
    + PartialEq
    + Eq
    + Default
    + ConditionallySelectable
{
    /// Create the constant zero
    fn zero() -> Self;

    /// Create the constant one  
    fn one() -> Self;

    /// Check if element is zero
    fn is_zero(&self) -> Choice;

    /// Check if element is one
    fn is_one(&self) -> Choice;

    /// Negate `self` to `-self`
    fn negate(&self) -> Self;

    /// Add `rhs` to `self`
    fn add(&self, rhs: &Self) -> Self;

    /// Sub `rhs` from `self`
    fn sub(&self, rhs: &Self) -> Self;

    /// Multipliy `self` by `rhs`
    fn mul(&self, rhs: &Self) -> Self;

    /// Square `self`
    fn square(&self) -> Self;

    /// Double `self`
    fn double(&self) -> Self;

    /// Invert `self`
    fn invert(&self) -> CtOption<Self>;

    /// Divide `self` by `rhs`
    fn div(&self, rhs: &Self) -> CtOption<Self> {
        rhs.invert().map(|inv| self.mul(&inv))
    }

    /// `self^exp` using square-and multiply (litte-endian bit order)
    ///
    /// # Arguments
    ///
    /// * `&self` - Finite field element (type: `Self`)
    /// * `exp` - Exponent (type: &[u64])
    ///
    /// # Returns
    ///
    /// `&self^exp` (type: `Self`)
    ///
    /// # Warnings
    ///
    /// This function is constant time only if two exponents `exp1`
    /// and `exp2` are equal as [u64] (i.e., they are the same number
    /// AND the same number of limbs)
    //
    // # Notes
    //
    // `<Self as FieldOps>::mul` is used instead of
    // `result.mul(&base)` because `FieldOps` requires `Mul<Output =
    // Self>` as a supertrait, so `Self` exposes **two** methods
    // named `mul`:
    //
    //   - `<Self as Mul>::mul(self, rhs: Self) -> Self` ← operator,
    //   takes by value - `<Self as FieldOps>::mul(&self, rhs: &Self)
    //   -> Self` ← ours, takes by ref
    //
    // Writing `result.mul(&base)` triggers method resolution, which
    // picks `Mul::mul` (the operator) because it was declared first
    // in the supertrait list. `Mul::mul` expects `Self`, not `&Self`
    // → E0308.
    //
    // Fully-qualified syntax `<Self as FieldOps>::mul(...)` bypasses
    // method resolution entirely and calls exactly the trait method
    // we want.
    fn pow_vartime(&self, exp: &[u64]) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();

        for &limb in exp {
            let mut limb = limb;
            for _ in 0..64 {
                let mybit = limb & 1;
                let tmp = <Self as FieldOps>::mul(&result, &base);
                result = Self::conditional_select(&result, &tmp, (mybit as u8).into());
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
    /// * `&self` - Element of $\mathbb{F}\_{p^M}$ (type: `Self`)
    /// * `exp` - Exponent (type: &[u64])
    ///
    /// # Returns
    ///
    /// The value `self^pow` (type: `Self`)
    ///
    /// # Warning
    ///
    /// This function is only constant time if the number of limbs of
    /// `exp` is fixed (so `[1,0]` would run in different time to
    /// `[1]`)
    fn pow(&self, exp: &[u64]) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();

        for &limb in exp.iter().rev() {
            let mut limb = limb.reverse_bits();
            for _ in 0..64 {
                let mybit = limb & 1;
                Self::conditional_swap(&mut base, &mut result, (mybit as u8).into());
                base = <Self as FieldOps>::mul(&result, &base);
                result = <Self as FieldOps>::square(&result);
                Self::conditional_swap(&mut base, &mut result, (mybit as u8).into());
                limb >>= 1;
            }
        }
        result
    }

    /// Compute `self^p` the frobenius acting on `self`
    fn frobenius(&self) -> Self;

    /// Compute `self^{p^k}` a power of the frobenius
    fn frobenius_pow(&self, k: u32) -> Self {
        let mut result = self.clone();
        for _ in 0..k {
            result = result.frobenius();
        }
        result
    }

    /// Compute the norm of `self` down to $\mathbb{F}_p$ (as an
    /// element of type `Self`)
    fn norm(&self) -> Self;

    /// Compute the trace of `self` down to $\mathbb{F}_p$ (as an
    /// element of type `Self`)
    fn trace(&self) -> Self;

    /// Returns a squareroot if it exists
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_{p^M}$ (type: `Self`)
    ///
    /// # Returns
    ///
    /// * $\sqrt{\texttt{self}}$ which is not `none` if and only if
    ///   the square root of `self` exists (type: `CtOption<Self>`)
    fn sqrt(&self) -> CtOption<Self>;

    /// Computes the inverse and square root of `self`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_{p^M}$ (type: `Self`)
    ///
    /// # Returns
    ///
    /// * `(inv, sqrt)` the inverse and the square root of `self`. The
    ///   former is none if and only if nonzero and the latter is not
    ///   none if and only if there exists a squareroot in
    ///   $\mathbb{F}_{p^M}$ (type: `(CtOption<Self>, CtOption<self>)`)
    fn inverse_and_sqrt(&self) -> (CtOption<Self>, CtOption<Self>) {
        (self.invert(), self.sqrt())
    }

    /// Computes the square root the inverse of `self`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_{p^M}$ (type: `Self`)
    ///
    /// # Returns
    ///
    /// * $1 / \sqrt{\texttt{self}}$. The is not none if and only if
    ///   `self` is both nonzero there exists a squareroot in
    ///   $\mathbb{F}_{p^M}$ (type: `CtOption<self>`)
    fn inv_sqrt(&self) -> CtOption<Self> {
        self.sqrt().and_then(|s| s.invert())
    }

    /// Computes the inverse of `self` and square root of `rhs`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_{p^M}$ (type: `Self`)
    /// * `rhs` - Element of $\mathbb{F}_{p^M}$ (type: &Self)
    ///
    /// # Returns
    ///
    /// * `(inv, sqrt)` where `inv` is the inverse of `self` `sqrt` is
    ///   the square root fo `rhs`. `inv` is `none` if and only if
    ///   `self` is zero and the `sqrt` is not none if and only if
    ///   there exists a squareroot of `rhs` in $\mathbb{F}_{p^M}$
    ///   (type: `(CtOption<Self>, CtOption<self>)`)
    fn invertme_sqrtother(&self, rhs: &Self) -> (CtOption<Self>, CtOption<Self>) {
        (self.invert(), rhs.sqrt())
    }

    /// Computes the squareroot of a ratio `self/rhs`
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of $\mathbb{F}_{p^M}$ (type: `Self`)
    /// * `rhs` - Element of $\mathbb{F}_{p^M}$ (type: `&Self`)
    ///
    /// # Returns
    ///
    /// * $\sqrt{\texttt{self} / \texttt{rhs}}$ which is not `none` if and only
    ///   both `rhs` is invertible and the ratio has a rational squareroot
    ///   (type: `(CtOption<Self>, CtOption<self>))`
    fn sqrt_ratio(&self, rhs: &Self) -> CtOption<Self> {
        rhs.invert().and_then(|inv_rhs| self.mul(&inv_rhs).sqrt())
    }

    /// Computes the "Legendre symbol" i.e., if 0,1,-1 depending if
    /// `self` is 0, a square or a nonsquare.
    fn legendre(&self) -> i8;

    /// Returns the characteristic of the field.
    fn characteristic() -> Vec<u64>;

    /// Returns the extension degree of the field.
    fn degree() -> u32;

    /// Convert u64 to the field.
    fn from_u64(x: u64) -> Self;
}
/// Trait for field types that can be constructed from a field-specific
/// representation.
pub trait FieldFromRepr: FieldOps {
    /// The representation type accepted by this field.
    type Repr;

    /// Constructs a field element from the given representation.
    fn from_repr(x: Self::Repr) -> Self;
}


// -------------------------------------------------------------------
// Macros for using operator overloads on referenced values
// -------------------------------------------------------------------

#[macro_export]
macro_rules! ref_field_impl {
    (impl<$F:ident> $Ty:ty { $($body:tt)* }) => {
        impl<$F> $Ty
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        {
            $($body)*
        }
    };

    (impl<$F:ident : $First:ident $(+ $Rest:ident)*> $Ty:ty { $($body:tt)* }) => {
        impl<$F: $First $(+ $Rest)*> $Ty
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        {
            $($body)*
        }
    };
}

#[macro_export]
macro_rules! ref_field_trait_impl {
    (impl<$F:ident> $Trait:ident for $Ty:ty { $($body:tt)* }) => {
        impl<$F> $Trait for $Ty
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        {
            $($body)*
        }
    };

    (impl<$F:ident : $First:ident $(+ $Rest:ident)*> $Trait:ident for $Ty:ty { $($body:tt)* }) => {
        impl<$F: $First $(+ $Rest)*> $Trait for $Ty
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        {
            $($body)*
        }
    };
}

#[macro_export]
macro_rules! ref_field_trait_impl_path {
    (impl<$F:ident> ($Trait:path) for $Ty:ty { $($body:tt)* }) => {
        impl<$F> $Trait for $Ty
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        {
            $($body)*
        }
    };

    (impl<$F:ident : $First:ident $(+ $Rest:ident)*> ($Trait:path) for $Ty:ty { $($body:tt)* }) => {
        impl<$F: $First $(+ $Rest)*> $Trait for $Ty
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        {
            $($body)*
        }
    };
}

#[macro_export]
macro_rules! ref_field_fn {
    (fn $name:ident<$F:ident> ( $($args:tt)* ) -> $Ret:ty $body:block) => {
        fn $name<$F>($($args)*) -> $Ret
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        $body
    };

    (fn $name:ident<$F:ident : $First:ident $(+ $Rest:ident)*> ( $($args:tt)* ) -> $Ret:ty $body:block) => {
        fn $name<$F: $First $(+ $Rest)*>($($args)*) -> $Ret
        where
            $F: $crate::field_ops::FieldOps,
            for<'a, 'b> &'a $F: std::ops::Add<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Sub<&'b $F, Output = $F>,
            for<'a, 'b> &'a $F: std::ops::Mul<&'b $F, Output = $F>,
            for<'a> &'a $F: std::ops::Neg<Output = $F>,
        $body
    };
}

#[macro_export]
macro_rules! ref_field_fns {
    () => {};

    (
        $(#[$meta:meta])*
        fn $name:ident<$F:ident> ( $($args:tt)* ) -> $Ret:ty $body:block
        $($rest:tt)*
    ) => {
        $(#[$meta])*
        $crate::ref_field_fn! {
            fn $name<$F>($($args)*) -> $Ret $body
        }
        $crate::ref_field_fns! { $($rest)* }
    };

    (
        $(#[$meta:meta])*
        fn $name:ident<$F:ident : $First:ident $(+ $Rest:ident)*> ( $($args:tt)* ) -> $Ret:ty $body:block
        $($rest:tt)*
    ) => {
        $(#[$meta])*
        $crate::ref_field_fn! {
            fn $name<$F: $First $(+ $Rest)*>($($args)*) -> $Ret $body
        }
        $crate::ref_field_fns! { $($rest)* }
    };
}