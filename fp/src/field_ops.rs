//! The [`FieldOps`](self::field_ops::FieldOps) trait is the core
//! trait that every field must implement.

use subtle::{Choice, ConditionallySelectable, CtOption};

/// Trait for generating cryptographically secure random field elements.
///
/// Separated from [`FieldOps`] so that downstream code that doesn't
/// need randomness is free of the `rand` dependency.
///
/// # Required Methods
///
/// * `random` - Generate a random field element.
///
/// # Examples
///
/// ```
/// use crypto_bigint::{const_prime_monty_params, Uint};
/// use fp::fp_element::FpElement;
/// use fp::field_ops::{FieldOps, FieldRandom};
///
/// /*
/// We will set up the field F_19 using FpElement, see the docs
/// of fp_element::FpElement.
/// */
///
/// const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2);
/// type F19 = FpElement<Fp19Modulus, 1>;
///
/// let mut rng = rand::rng();
/// let a = F19::random(&mut rng);
/// ```
pub trait FieldRandom: Sized {
    /// Sample a uniformly random element using a cryptographic random
    /// number generator.
    ///
    /// # Arguments
    ///
    /// * `rng` a cryptographic rng (type: `&mut (impl rand::CryptoRng + rand::Rng)`)
    ///
    /// # Returns
    ///
    /// A random element of $\mathbb{F}\_{q}$ (type: `Self`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::{FieldOps, FieldRandom};
    /// let mut rng = rand::rng();
    /// let a = F19::random(&mut rng);
    /// ```
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self;
}

/// Core arithmetic interface for field elements used throughout the library.
///
/// This trait abstracts the operations needed by the prime-field layer,
/// extension fields, and higher-level elliptic-curve code.
///
/// # Required Methods
///
/// * `zero` - Returns the constant zero.
/// * `one` - Returns the constant one.
/// * `from_u64` - Coerces as `u64` into the field.
/// * `is_zero` - Checks if element is zero using constant time `Choice`.
/// * `is_one` - Checks if element is zero using constant time `Choice`.
/// * `negate` - Returns $-\texttt{self}$.
/// * `add` - Returns $\texttt{self} + \texttt{rhs}$.
/// * `sub` - Returns $\texttt{self} - \texttt{rhs}$.
/// * `mul` - Returns $\texttt{self} * \texttt{rhs}$.
/// * `square` - Returns $\texttt{self}^2$.
/// * `double` - Returns $2 * \texttt{self}$.
/// * `invert` - Returns $\texttt{self}^{-1}$.
/// * `frobenius` - Returns $\texttt{self}^{\mathrm{char} \mathbb{F}\_{q}}$.
/// * `norm` - Returns $\mathrm{Norm}\_{\mathbb{F}\_{q} /
///   \mathbb{F}\_{p}} \texttt{self}$.
/// * `trace` - Returns $\mathrm{Tr}\_{\mathbb{F}\_{q} /
///   \mathbb{F}\_{p}} \texttt{self}$.
/// * `sqrt` - Returns $\sqrt{\texttt{self}}$.
/// * `legendre` - Returns the "Legendre symbol" which is 1 if and
///   only if $\texttt{self}$ is square.
/// * `characteristic` - Returns the characteristic of the field.
/// * `degree` - Returns $[\mathbb{F}\_{p^M} : \mathbb{F}_p]$
///
/// # Provided Methods
///
/// * `div` - Returns $\texttt{self} / \texttt{rhs}$.
/// * `pow_vartime` - Returns $\texttt{self}^\texttt{pow}$ using square and
///   multiply.
/// * `pow` - Computes $\texttt{self}^\texttt{pow}$ in constant time using a
///   Montgomery ladder.
/// * `frobenius_pow` - Computes $\texttt{self}^{p^k}$
/// * `inverse_and_sqrt` - Computes $\texttt{self}^{-1}$ and $\sqrt{\texttt{self}}$.
/// * `inv_sqrt` - Computes $\sqrt{\texttt{self}^{-1}}$.
/// * `invertme_sqrtother` - Computes $\texttt{self}^{-1}$ and $\sqrt{\texttt{rhs}}$.
/// * `sqrt_ratio` - Computes $\sqrt{\texttt{self} / \texttt{rhs}}$.
///
/// Scalars used by `pow_vartime` and `pow` are encoded as
/// little-endian `u64` limbs, matching the convention used elsewhere
/// in the library.
///
/// # Note
///
/// The trait is intentionally separate from [`FieldRandom`], so downstream code
/// that only needs deterministic arithmetic does not need to depend on `rand`.
///
/// # Examples
///
/// ## The field $\mathbb{F}\_{19}$
///
/// ```
/// use crypto_bigint::{const_prime_monty_params, Uint};
/// use fp::fp_element::FpElement;
/// use fp::field_ops::FieldOps;
///
/// /*
/// We will set up the field F_19 using FpElement, see the docs
/// of fp_element::FpElement.
/// */
///
/// const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2);
/// type F19 = FpElement<Fp19Modulus, 1>;
///
/// /* Some standard tests */
/// let a = F19::from_u64(17);
/// let b = F19::from_u64(5);
/// assert_eq!(a.add(&b), F19::from_u64(3));
/// assert_eq!(a.sub(&b), F19::from_u64(12));
/// assert_eq!(a.mul(&b), F19::from_u64(9));
/// assert_eq!(a.invert().unwrap(), F19::from_u64(9));
/// ```
///
/// ## The field $\mathbb{F}\_{19^2}$
///
/// ```
/// # use crypto_bigint::{const_prime_monty_params, Uint};
/// # use fp::fp_element::FpElement;
/// # use fp::field_ops::FieldOps;
/// # const_prime_monty_params!(Fp19Modulus, Uint<1>, "0000000000000013", 2);
/// # type F19 = FpElement<Fp19Modulus, 1>;
/// /* Continuing from the above example */
/// use fp::fp_ext::{FpExt, IrreduciblePoly, TonelliShanksConstants};
///
/// /* Setput the irreducible polynomial */
/// struct QuadPoly;
/// impl IrreduciblePoly<Fp19Modulus, 1, 2> for QuadPoly {
///     fn modulus() -> [F19; 2] {
///         [F19::one(), F19::zero()]
///     }
/// }
///
/// /* Setup the Tonelli--Shanks constants */
/// struct TSQuad;
/// impl TonelliShanksConstants<Fp19Modulus, 1, 2, 1> for TSQuad {
///     // p^2 - 1
///     const ORDER: Uint<1> = Uint::<1>::from_u64(360);
///     // (p^2 - 1) / 2
///     const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(180);
///     // p^2 - 1 = 2^S * T with T odd
///     const S: u64 = 3;
///     const T: Uint<1> = Uint::<1>::from_u64(45);
///     // (T - 1) / 2
///     const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(22);
///     // 2^(S - 1)
///     const TWOSM1: Uint<1> = Uint::<1>::from_u64(4);
///     // 2^S root of unity
///     fn root_of_unity() -> [FpElement<Fp19Modulus, 1>; 2] {
///         [F19::from_u64(3), F19::from_u64(3)]
///     }
/// }
///
/// /* Define the extension field */
/// type F19_2 = FpExt<Fp19Modulus, 1, 2, 1, QuadPoly, TSQuad>;
///
/// /* Some standard tests */
/// /* Some standard tests */
/// let a = F19_2::from_u64(17);
/// let b = F19_2::from_u64(5);
/// assert_eq!(a.add(&b), F19_2::from_u64(3));
/// assert_eq!(a.sub(&b), F19_2::from_u64(12));
/// assert_eq!(a.mul(&b), F19_2::from_u64(9));
/// assert_eq!(a.invert().unwrap(), F19_2::from_u64(9));
/// ```
pub trait FieldOps: Sized + Clone + PartialEq + Eq + Default + ConditionallySelectable {
    /// Create the constant zero
    ///
    /// # Returns
    ///
    /// The constant $0 \in \mathbb{F}\_q$ (type: `Self`).
    fn zero() -> Self;

    /// Create the constant one
    ///
    /// # Returns
    ///
    /// The constant $1 \in \mathbb{F}\_q$ (type: `Self`).
    fn one() -> Self;

    /// Convert `u64` to the field.
    ///
    /// # Arguments
    ///
    /// * `x` an 64 bit integer (type: `u64`).
    ///
    /// # Returns
    ///
    /// The element $x \in \mathbb{F}\_p$ embedded in $\mathbb{F}\_q$
    /// (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let my_zero = F19::from_u64(0);
    /// assert_eq!(my_zero, F19::zero());
    /// let my_one = F19::from_u64(1);
    /// assert_eq!(my_one, F19::one());    
    /// ```
    fn from_u64(x: u64) -> Self;

    /// Checks if an element is zero
    ///
    /// # Returns
    ///
    /// `true` if and only if `self` is zero (type: `Choice`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let x = F19::zero();
    /// assert!(bool::from(x.is_zero()));
    /// ```
    fn is_zero(&self) -> Choice;

    /// Checks if an element is one
    ///
    /// # Returns
    ///
    /// `true` if and only if `self` is one (type: `Choice`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let x = F19::one();
    /// assert!(bool::from(x.is_one()));
    /// ```
    fn is_one(&self) -> Choice;

    /// Negate `self` to `-self`
    ///
    /// # Returns
    ///
    /// $-\texttt{self}$ (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let x = F19::from_u64(17);
    /// assert_eq!(x.negate(), F19::from_u64(2));
    /// ```
    fn negate(&self) -> Self;

    /// Returns $\texttt{self} + \texttt{rhs}$.
    ///
    /// # Arguments
    ///
    /// * `rhs` element of $\mathbb{F}\_{q}$ (type: `&Self`).
    ///
    /// # Returns
    ///
    /// $\texttt{self} + \texttt{rhs}$ (type: `Self`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(10);
    /// let b = F19::from_u64(11);
    /// let a_add_b = a.add(&b);
    /// assert_eq!(a_add_b, F19::from_u64(2));
    /// ```
    fn add(&self, rhs: &Self) -> Self;

    /// Returns $\texttt{self} - \texttt{rhs}$.
    ///
    /// # Arguments
    ///
    /// * `rhs` element of $\mathbb{F}\_{q}$ (type: `&Self`).
    ///
    /// # Returns
    ///
    /// $\texttt{self} - \texttt{rhs}$ (type: `Self`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(4);
    /// let b = F19::from_u64(8);
    /// let a_sub_b = a.sub(&b);
    /// assert_eq!(a_sub_b, F19::from_u64(15));
    /// ```
    fn sub(&self, rhs: &Self) -> Self;

    /// Returns $\texttt{self} * \texttt{rhs}$.
    ///
    /// # Arguments
    ///
    /// * `rhs` element of $\mathbb{F}\_{q}$ (type: `&Self`).
    ///
    /// # Returns
    ///
    /// $\texttt{self} * \texttt{rhs}$ (type: `Self`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(4);
    /// let b = F19::from_u64(3);
    /// let a_times_b = a.mul(&b);
    /// assert_eq!(a_times_b, F19::from_u64(12));
    /// ```
    fn mul(&self, rhs: &Self) -> Self;

    /// Returns $\texttt{self}^2$.
    ///
    /// # Returns
    ///
    /// $\texttt{self}^2$ (type: `Self`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(4);
    /// let a_sq = a.square();
    /// assert_eq!(a_sq, F19::from_u64(16));
    /// ```
    fn square(&self) -> Self;

    /// Returns $2 * \texttt{self}$.
    ///
    /// # Returns
    ///
    /// $2 * \texttt{self}$ (type: `Self`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(4);
    /// let a_doub = a.double();
    /// assert_eq!(a_doub, F19::from_u64(8));
    /// ```
    fn double(&self) -> Self;

    /// Returns $\texttt{self}^{-1}$.
    ///
    /// # Returns
    ///
    /// $\texttt{self}^{-1}$ which is none if and only if `self` is
    /// nonzero (type: `CtOption<Self>`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(2);
    /// let a_inv = a.invert().unwrap();
    /// assert_eq!(a_inv, F19::from_u64(10));
    /// let z = F19::zero();
    /// let z_inv = z.invert();
    /// assert!(bool::from(z_inv.is_none()));
    /// ```
    fn invert(&self) -> CtOption<Self>;

    /// Returns $\texttt{self} / \texttt{rhs}$.
    ///
    /// # Arguments
    ///
    /// * `rhs` element of $\mathbb{F}\_{q}$ (type: `&Self`).
    ///
    /// # Returns
    ///
    /// $\texttt{self} / \texttt{rhs}$ which is none if and only if
    /// `rhs` is nonzero (type: `CtOption<Self>`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(4);
    /// let b = F19::from_u64(2);
    /// let a_div_b = a.div(&b).unwrap();
    /// assert_eq!(a_div_b, F19::from_u64(2));
    /// let z = F19::zero();
    /// let fail = a.div(&z);
    /// assert!(bool::from(fail.is_none()));
    /// ```
    fn div(&self, rhs: &Self) -> CtOption<Self> {
        rhs.invert().map(|inv| self.mul(&inv))
    }

    /// Computes $\texttt{self}^\texttt{exp}$ in vartime.
    ///
    /// # Arguments
    ///
    /// * `exp` - Exponent (type: `&[u64]`)
    ///
    /// # Returns
    ///
    /// `self^exp` (type: `Self`)
    ///
    /// # WARNING
    ///
    /// This function is constant time only if two exponents `exp1`
    /// and `exp2` are equal as [u64] (i.e., they are the same number
    /// AND the same number of limbs)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(3);
    /// let a_cubed = a.pow_vartime(&[3]);
    /// assert_eq!(a_cubed, F19::from_u64(8));
    /// ```
    ///
    /// # Notes
    ///
    /// The provided function is implemented using square-and
    /// multiply (litte-endian bit order).
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

    /// Computes $\texttt{self}^\texttt{exp}$ in constant time.
    ///
    /// # Arguments
    ///
    /// * `exp` - Exponent (type: `&[u64]`)
    ///
    /// # Returns
    ///
    /// `self^exp` (type: `Self`)
    ///
    /// # WARNING
    ///
    /// This function is only required to be constant time if the
    /// number of limbs of `exp` is fixed (so `[1,0]` would run in
    /// different time to `[1]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(3);
    /// let a_cubed = a.pow(&[3]);
    /// assert_eq!(a_cubed, F19::from_u64(8));
    /// ```
    ///
    /// # Note
    ///
    /// The provided function is implemented using a Montgomery ladder
    /// and `mul`.
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

    /// Compute $\texttt{self}^p$ the Frobenius endomorphism
    ///
    /// # Returns
    ///
    /// $\texttt{self}^p$ (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// assert_eq!(a.frobenius(), a);
    /// let b = F19_2::from_u64_vec([0,1]); // The element i with i^2 = -1
    /// assert_eq!(b.frobenius(), b.negate());
    /// ```
    fn frobenius(&self) -> Self;

    /// Compute $\texttt{self}^{p^k}$ a power of the Frobenius
    /// endomorphism
    ///
    /// # Arguments
    ///
    /// * `k` a positive integer (type: `u32`).
    ///
    /// # Returns
    ///
    /// $\texttt{self}^{p^k}$ (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// assert_eq!(a.frobenius_pow(3), a);
    /// let b = F19_2::from_u64_vec([0,1]); // The element i with i^2 = -1
    /// assert_eq!(b.frobenius_pow(2), b);
    /// assert_eq!(b.frobenius_pow(3), b.frobenius());
    /// ```
    fn frobenius_pow(&self, k: u32) -> Self {
        let mut result = self.clone();
        for _ in 0..k {
            result = result.frobenius();
        }
        result
    }

    /// Computes the norm of `self` to $\mathbb{F}\_{p}$ embedded as an element
    /// of $\mathbb{F}\_{q}$.
    ///
    /// As a reminder the norm is defined as the product over all the
    /// Galois conjugates that is
    /// $$ N\_{\mathbb{F}\_{q}/\mathbb{F}\_p}(a) = \prod\_{\sigma \in \mathrm{Gal}} \sigma(a) $$
    ///
    /// # Returns
    ///
    /// The norm of self (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19_2::from_u64_vec([3, 5]);
    /// let nm_a = F19_2::from_u64_vec([15, 0]);
    /// assert_eq!(a.norm(), nm_a);
    /// ```
    fn norm(&self) -> Self;

    /// Computes the trace of `self` to $\mathbb{F}\_{p}$ embedded as an element
    /// of $\mathbb{F}\_{q}$.
    ///
    /// As a reminder the trace is defined as the sum over all the
    /// Galois conjugates that is
    /// $$ \mathrm{tr}\_{\mathbb{F}\_{q}/\mathbb{F}\_p}(a) = \sum\_{\sigma \in \mathrm{Gal}} \sigma(a) $$
    ///
    /// # Returns
    ///
    /// The trace of self (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19_2::from_u64_vec([3, 5]);
    /// let tr_a = F19_2::from_u64_vec([6, 0]);
    /// assert_eq!(a.trace(), tr_a);
    /// ```
    fn trace(&self) -> Self;

    /// Computes the square root of `self`
    ///
    /// # Returns
    ///
    /// $\sqrt{\texttt{self}}$ a choice of squareroot which is not
    /// `none` if the squareroot exists (type: `CtOption<Self>`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// let sqrt_a = a.sqrt().unwrap();
    /// assert_eq!(sqrt_a.mul(&sqrt_a), a);
    /// ```
    fn sqrt(&self) -> CtOption<Self>;

    /// Returns both the inverse and sqrt of `self`
    ///
    /// # Returns
    ///
    /// $(1/\texttt{self}, \sqrt{\texttt{self}})$ which is
    /// `self.invert()` and `self.sqrt()` (type: `CtOption<Self>`,
    /// `CtOption<Self>`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// let (inv_a, sqrt_a) = a.inverse_and_sqrt();
    /// let inv_a = inv_a.unwrap();
    /// let sqrt_a = sqrt_a.unwrap();
    /// assert!(bool::from(inv_a.mul(&a).is_one()));
    /// assert_eq!(sqrt_a.mul(&sqrt_a), a);
    /// ```
    ///
    /// # Note
    ///
    /// The default implementation just combines the `sqrt` and
    /// `invert` methods. This is provided as a separate method so
    /// that a motivated user can implement tricks such as those found
    /// in Section 2 of <https://eprint.iacr.org/2020/1497>. See e.g.,
    /// implementation in [`FpExt`](crate::fp_ext::FpExt).
    fn inverse_and_sqrt(&self) -> (CtOption<Self>, CtOption<Self>) {
        (self.invert(), self.sqrt())
    }

    /// Computes the inverse and squareroot of `self`.
    ///
    /// # Returns
    ///
    /// $1 / \sqrt{\texttt{self}}$, the inverse of the squareroot of
    /// `self` (type: `CtOption<Self>`)
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// let inv_sqrt_a = a.inv_sqrt();
    /// let inv_sqrt_a = inv_sqrt_a.unwrap();
    /// assert_eq!(inv_sqrt_a.mul(&inv_sqrt_a), a.invert().unwrap());
    /// ```
    ///
    /// # Note
    ///
    /// The default implementation just combines the `sqrt` and
    /// `invert` methods. This is provided as a separate method so
    /// that a motivated user can implement tricks such as those found
    /// in Section 2 of <https://eprint.iacr.org/2020/1497>. See e.g.,
    /// implementation in [`FpExt`](crate::fp_ext::FpExt).
    fn inv_sqrt(&self) -> CtOption<Self> {
        self.sqrt().and_then(|s| s.invert())
    }

    /// Inverse of `self` and squareroot of `rhs`.
    ///
    /// # Arguments
    ///
    /// * `rhs` - element of $\mathbb{F}\_{p^M}$ (type: `&Self`)
    ///
    /// # Returns
    ///
    /// $(1 / \texttt{self}, \sqrt{\texttt{rhs}})$, the former is none
    /// if and only if `self` is nonzero and the latter is none if and
    /// only if there is no squareroot of `rhs` in $\mathbb{F}\_{p^M}$
    /// (type: (`CtOption<Self>`, `CtOption<Self>`))
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// let b = F19::from_u64(4);
    /// let (inv_a, sqrt_b) = a.invertme_sqrtother(&b);
    /// let inv_a = inv_a.unwrap();
    /// let sqrt_b = sqrt_b.unwrap();
    /// assert!(bool::from(inv_a.mul(&a).is_one()));
    /// assert_eq!(sqrt_b.square(), b);
    /// ```
    ///
    /// # Note
    ///
    /// The default implementation just combines the `sqrt` and
    /// `invert` methods. This is provided as a separate method so
    /// that a motivated user can implement tricks such as those found
    /// in Section 2 of <https://eprint.iacr.org/2020/1497>. See e.g.,
    /// implementation in [`FpExt`](crate::fp_ext::FpExt).
    fn invertme_sqrtother(&self, rhs: &Self) -> (CtOption<Self>, CtOption<Self>) {
        (self.invert(), rhs.sqrt())
    }

    /// Computes the squareroot of a ratio `self/rhs`
    ///
    /// # Arguments
    /// * `rhs` - Element of $\mathbb{F}\_{p^M}$ (type: `&Self`)
    ///
    /// # Returns
    ///
    /// $\sqrt{\texttt{self}/\texttt{rhs}}$, the squareroot of the
    /// ratio which is not none if and only if `rhs` is invertible and
    /// the ratio has an $\mathbb{F}\_{p^M}$ rational squareroot (type:
    /// `CtOption<Self>`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// let b = F19::from_u64(4);
    /// let sqrt_a_over_b = a.sqrt_ratio(&b);
    /// let sqrt_a_over_b = sqrt_a_over_b.unwrap();
    /// assert_eq!(sqrt_a_over_b.square().mul(&b), a);
    /// ```
    ///
    /// # Note
    ///
    /// The default implementation just combines the `sqrt` and
    /// `invert` methods. This is provided as a separate method so
    /// that a motivated user can implement tricks such as those found
    /// in Section 2 of <https://eprint.iacr.org/2020/1497>. See e.g.,
    /// implementation in [`FpExt`](crate::fp_ext::FpExt).    
    fn sqrt_ratio(&self, rhs: &Self) -> CtOption<Self> {
        rhs.invert().and_then(|inv_rhs| self.mul(&inv_rhs).sqrt())
    }

    /// Implements the "Legendre symbol" which is 1 if and only if
    /// `self` is a quadratic residue in $\mathbb{F}\_{q}$.
    ///
    /// As a reminder the "Legendre symbol" for $a \in \mathbb{F}\_{q}$ is defined as
    /// $$ \begin{cases} 0 & \text{ if $a = 0$,} \\\\ 1 & \text{if $a$
    /// is a square,} \\\\ -1 & \text{if $a$ is not a square.} \end{cases} $$
    ///
    /// # Returns
    ///
    /// The Legendre symbol of `self` (type: `i8`)
    ///
    /// # WARNING
    ///
    /// Not required to be constant time if `self` is zero
    ///
    /// # Examples
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let a = F19::from_u64(9);
    /// let b = F19::from_u64(2);
    /// let c = F19::zero();
    /// let ls_a = a.legendre();
    /// let ls_b = b.legendre();
    /// let ls_c = c.legendre();
    /// assert_eq!(ls_a, 1_i8);
    /// assert_eq!(ls_b, -1_i8);
    /// assert_eq!(ls_c, 0_i8);
    /// ```
    fn legendre(&self) -> i8;

    /// Returns the characteristic of the field.
    ///
    /// # Returns
    ///
    /// The characteristic of $\mathbb{F}\_{q}$ as 64 bit limbs (type:
    /// `Vec<u64>`)
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let p = F19::characteristic();
    /// let p_2 = F19_2::characteristic();
    /// assert_eq!(p, p_2);
    /// assert_eq!(p[0], 19_u64);
    /// ```
    fn characteristic() -> Vec<u64>;

    /// Returns the degree of the field over $\mathbb{F}_p$.
    ///
    /// # Returns
    ///
    /// The degree of $\mathbb{F}\_{q} / \mathbb{F}\_{p}$ (type:
    /// `u32`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_field_ops::*;
    /// # use fp::field_ops::FieldOps;
    /// let d_1 = F19::degree();
    /// let d_2 = F19_2::degree();
    /// assert_eq!(d_1, 1_u32);
    /// assert_eq!(d_2, 2_u32);
    /// ```
    fn degree() -> u32;
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

/// Helper macro for generic code that wants to use borrowed operators
/// like `&a + &b`, `&a - &b`, `&a * &b`, and `-&a`.
///
/// # Example
///
/// ```rust
/// use fp::field_ops::FieldOps;
/// use fp::ref_field_impl;
///
/// pub struct PairOfFieldElts<F: FieldOps> {
///     a: F,
///     b: F,
/// }
///
/// ref_field_impl! {
///     impl<F> PairOfFieldElts<F> {
///         pub fn new(a: F, b: F) -> Self {
///             Self { a, b }
///         }
///
///         pub fn add(&self) -> F {
///             &self.a + &self.b
///         }
///     }
/// }
/// ```
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

/// Helper macro for trait implementations that wants to use borrowed
/// operators like `&a + &b`, `&a - &b`, `&a * &b`, and `-&a`.
///
/// # Example
///
/// ```
/// use fp::field_ops::FieldOps;
/// use fp::ref_field_trait_impl;
///
/// pub struct PairOfFieldElts<F: FieldOps> {
///     a: F,
///     b: F,
/// }
///
/// pub trait Data {
///     type BaseField;
///     fn are_equal(&self) -> bool;
/// }
///
/// ref_field_trait_impl! {
///     impl<F: FieldOps> Data for PairOfFieldElts<F> {
///         type BaseField = F;
///         fn are_equal(&self) -> bool {
///             self.a == self.b
///         }
///     }
/// }
/// ```
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

/// For trait impls when one wants to specify the path and also use
/// borrowed operators like `&a + &b`, `&a - &b`, `&a * &b`, and
/// `-&a`.
///
/// # Example
///
/// ```
/// use fp::field_ops::FieldOps;
/// use fp::ref_field_trait_impl_path;
///
/// pub struct PairOfFieldElts<F: FieldOps> {
///     a: F,
///     b: F,
/// }
///
/// # mod path {
/// #     pub mod to {
/// #         pub mod dependency {
/// #             pub trait Data {
/// #                 type BaseField;
/// #                 fn are_equal(&self) -> bool;
/// #             }
/// #         }
/// #     }
/// # }
/// #
/// ref_field_trait_impl_path! {
///     impl<F> (path::to::dependency::Data) for PairOfFieldElts<F> {
///         type BaseField = F;
///         fn are_equal(&self) -> bool {
///             self.a == self.b
///         }
///     }
/// }
/// ```
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

/// Helper macro for functions which wish to use operator overloads
///
/// # Example
///
/// ```
/// use fp::field_ops::FieldOps;
/// use fp::ref_field_fn;
///
/// ref_field_fn! {
///     fn add_elts<F>(a1: &F, a2: &F) -> F {
///         a1 + a2
///     }
/// }
/// ```
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

/// For batches of functions which wish to use operator overloads
///
/// # Example
///
/// ```
/// use fp::field_ops::FieldOps;
/// use fp::ref_field_fns;
///
/// ref_field_fns! {
///     fn add_elts<F>(a1: &F, a2: &F) -> F {
///         a1 + a2
///     }
///
///     fn mul_elts<F>(a1: &F, a2: &F) -> F {
///         a1 * a2
///     }
/// }
/// ```
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
