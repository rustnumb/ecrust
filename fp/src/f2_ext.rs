//! Generic binary fields $\mathbb{F}\_{2^m} = \mathbb{F}\_2\[x\] / (f(x))$
//!
//! # Examples
//!
//! ```
//! use crypto_bigint::Uint;
//! use fp::f2_element::F2Element;
//! use fp::f2_ext::{BinaryIrreducible, F2Ext};
//! use fp::field_ops::FieldOps;
//!
//! /* Make the finite field F_4 */
//! struct F4Poly;
//! impl BinaryIrreducible<1> for F4Poly {
//!    fn modulus() -> Uint<1> {
//!        Uint::<1>::from_u64(0b111) // x^2 + x + 1
//!    }
//!
//!    fn degree() -> usize {
//!        2usize
//!    }
//! }
//! type F4 = F2Ext<1, F4Poly>;
//!
//! /* Elements are written down by binary integers */
//! let zero = F4::from_u64(0);
//! let one = F4::from_u64(1);
//! let a = F4::from_u64(0b11);
//! let b = F4::from_u64(0b10);
//!
//! let also_zero = F4::from_u64(0b111); // x^2 + x + 1 = 0
//! let also_one = F4::from_u64(0b1000); // x^3 = x*x^2 = x^2 + x = 1
//! assert_eq!(also_zero, zero);
//! assert_eq!(also_one, one);
//! assert_eq!(a.mul(&b), one);
//! ```

use crate::field_ops::{FieldFromRepr, FieldOps, FieldRandom};
use core::ops::{Add, Mul, Neg, Sub};
use crypto_bigint::Uint;
use std::marker::PhantomData;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// ---------------------------------------------------------------------------
// IrreduciblePoly — the only thing callers need to implement for a new field
// ---------------------------------------------------------------------------

/// Irreducible polynomial defined over $\mathbb{F}_{2}$
///
/// # Required Methods
///
/// * `modulus` - The polynomial defining the extension field
///   implemented as a `Uint<LIMBS>` one box for each coefficient
/// * `degree` - The degree of the extension.
///
/// # Examples
///
/// ```
/// use crypto_bigint::Uint;
/// use fp::f2_element::F2Element;
/// use fp::f2_ext::{BinaryIrreducible, F2Ext};
///
/// struct MyPoly;
///
/// impl BinaryIrreducible<1> for MyPoly {
///    fn modulus() -> Uint<1> {
///        Uint::<1>::from_u64(0b1_0001_1011) // x^8 + x^4 + x^3 + x + 1
///    }
///
///    fn degree() -> usize {
///        8usize
///    }
/// }
/// ```
pub trait BinaryIrreducible<const LIMBS: usize>: 'static {
    /// Full polynomial bitmask, including the leading term $x^m$
    ///
    /// # Returns
    ///
    /// The irreducible polynomial (type: `Uint<LIMBS>`).
    fn modulus() -> Uint<LIMBS>;

    /// Degree m of the irreducible polynomial
    ///
    /// # Returns
    ///
    /// The degree of the irreducible polynomial `modulus` (type:
    /// `usize`)
    fn degree() -> usize;
}

// ---------------------------------------------------------------------------
// F2Ext — element of F_{2^M}
// ---------------------------------------------------------------------------

/// An extension of $\mathbb{F}_2$ given by an irreducible binary
/// polynomial $P \in \mathbb{F}_2\[x\]$.
///
/// # Examples
///
/// ```
/// use crypto_bigint::Uint;
/// use fp::f2_element::F2Element;
/// use fp::f2_ext::{BinaryIrreducible, F2Ext};
///
/// /* Make the finite field F_4 */
/// struct F4Poly;
/// impl BinaryIrreducible<1> for F4Poly {
///    fn modulus() -> Uint<1> {
///        Uint::<1>::from_u64(0b111) // x^2 + x + 1
///    }
///
///    fn degree() -> usize {
///        2usize
///    }
/// }
/// type F4 = F2Ext<1, F4Poly>;
///
/// /* Make the finite field F_{2^511} */
/// struct F2511Poly;
/// impl BinaryIrreducible<8> for F2511Poly {
///    fn modulus() -> Uint<8> {
///        let one = Uint::<8>::from_u64(1);
///        (one << 511) | (one << 10) | one
///    }
///
///    fn degree() -> usize {
///        511
///    }
/// }
/// type F2_511 = F2Ext<8, F2511Poly>;
/// ```
pub struct F2Ext<const LIMBS: usize, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    /// The value of an element of $\mathbb{F}_{2^M}$
    pub value: Uint<LIMBS>,
    _phantom: PhantomData<P>,
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<const LIMBS: usize, P> F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    /// Make an element of the extension from the limbs
    ///
    /// # Arguments
    ///
    /// * `a` - An integer written in binary as $a_0 \dots a_{64 *
    ///   \texttt{LIMBS}}$ and padded with 0s if nessicary (type:
    ///   `Uint<LIMBS>`).
    ///
    /// # Returns
    ///
    /// The element $\sum_{i=0}^M a_i x^i \in \mathbb{F}\_{2}\[x\] /
    /// (f(x)) = \mathbb{F}\_{2^M}$ (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_f2_ext::*;
    /// # use crypto_bigint::Uint;
    /// # use fp::field_ops::FieldOps;
    /// let a = F4::new(Uint::<1>::from_u64(0b111)); // x^2 + x + 1 should be zero
    /// let b = F4::new(Uint::<1>::from_u64(0b11)); // x + 1 = x^2
    /// let c = F4::new(Uint::<1>::from_u64(0b100)); // x^2 = x + 1
    /// let d = F4::new(Uint::<1>::from_u64(0b1000)); // x^3 = x^2 * x = x^2 + x = 1
    /// assert!(bool::from(a.is_zero()));
    /// assert_eq!(b, c);
    /// assert!(bool::from(d.is_one()));
    /// ```
    pub fn new(a: Uint<LIMBS>) -> Self {
        Self {
            value: reduce::<LIMBS, P>(a),
            _phantom: PhantomData,
        }
    }

    /// Make an element of the extension from a `u64`
    ///
    /// # Arguments
    ///
    /// * `a` - An integer written in binary as $a_0 \dots a_{64}$ and
    /// padded with 0s to 64 bits if nessicary (type: `u64`).
    ///
    /// # Returns
    ///
    /// The element $\sum_{i=0}^{64} a_i x^i \in \mathbb{F}\_{2}\[x\] /
    /// (f(x)) = \mathbb{F}\_{2^M}$ (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_f2_ext::*;
    /// # use crypto_bigint::Uint;
    /// # use fp::field_ops::FieldOps;
    /// let a = F4::from_u64(0b111); // x^2 + x + 1 should be zero
    /// let b = F4::from_u64(0b11); // x + 1 should be x^2
    /// let c = F4::from_u64(0b100); // x + 1 should be x^2
    /// let d = F4::from_u64(0b1000); // x^3 = x^2 * x = x^2 + x = 1
    /// assert!(bool::from(a.is_zero()));
    /// assert_eq!(b, c);
    /// assert!(bool::from(d.is_one()));
    /// ```
    pub fn from_u64(x: u64) -> Self {
        Self::new(Uint::from(x))
    }

    /// Make an element of the extension from the limbs
    ///
    /// # Arguments
    ///
    /// * `a` - An integer written in binary as $a_0 a_1 a_2 \dots
    ///   a_{M}$ and then padded to the nearest 64 bits (type:
    ///   `Uint<LIMBS>`).
    ///
    /// # Returns
    ///
    /// The element $\sum_{i=0}^M a_i x^i \in \mathbb{F}\_{2}\[x\] /
    /// (f(x)) = \mathbb{F}\_{2^M}$ (type: `Self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_f2_ext::*;
    /// # use crypto_bigint::Uint;
    /// # use fp::field_ops::FieldOps;
    /// let a = F4::from_uint(Uint::<1>::from_u64(0b111)); // x^2 + x + 1 should be zero
    /// let b = F4::from_uint(Uint::<1>::from_u64(0b11)); // x + 1 = x^2
    /// let c = F4::from_uint(Uint::<1>::from_u64(0b100)); // x^2 = x + 1
    /// let d = F4::from_uint(Uint::<1>::from_u64(0b1000)); // x^3 = x^2 * x = x^2 + x = 1
    /// assert!(bool::from(a.is_zero()));
    /// assert_eq!(b, c);
    /// assert!(bool::from(d.is_one()));
    /// ```
    ///
    /// # Note
    ///
    /// This is just an alias for [`F2Ext::new`]
    pub fn from_uint(x: Uint<LIMBS>) -> Self {
        Self::new(x)
    }

    /// Get a `Uint<LIMBS>` from an element, inverting [`F2Ext::new`].
    ///
    /// # Returns
    ///
    /// The integer $a$ of at most $M$ bits such that applying `new`
    /// gives `self`. That is, if $a$ is written in binary as $a_0
    /// \dots a_M$ then $\texttt{self} = \sum_{i=0}^M a_i x^i \in
    /// \mathbb{F}\_{2}\[x\] / (f(x)) = \mathbb{F}\_{2^M}$ (type:
    /// `Uint<LIMBS>`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_f2_ext::*;
    /// # use crypto_bigint::Uint;
    /// # use fp::field_ops::FieldOps;
    /// let a = F4::from_uint(Uint::<1>::from_u64(0b11)); // x + 1
    /// let b = F4::from_uint(Uint::<1>::from_u64(0b1000)); // x^3 = x^2 * x = x^2 + x = 1
    /// assert_eq!(a.as_uint(), Uint::<1>::from_u64(0b11));
    /// assert_eq!(b.as_uint(), Uint::<1>::from_u64(0b1));
    /// ```
    pub fn as_uint(&self) -> Uint<LIMBS> {
        self.value
    }

    /// Get the degree of the field extension
    ///
    /// # Returns
    ///
    /// The degree of the field extension (type: `usize`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use fp::_doctest_support::_doctest_f2_ext::*;
    /// # use crypto_bigint::Uint;
    /// # use fp::field_ops::FieldOps;
    /// let d = F4::degree();
    /// assert_eq!(d, 2_usize);
    /// ```
    pub fn degree() -> usize {
        P::degree()
    }
}

// ---------------------------------------------------------------------------
// Clone / Copy / PartialEq / Eq / Debug
// (manual impls so we don't over-constrain the bounds the way #[derive] would)
// ---------------------------------------------------------------------------

impl<const LIMBS: usize, P> Clone for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn clone(&self) -> Self {
        Self {
            value: self.value.clone(),
            _phantom: PhantomData,
        }
    }
}

impl<const LIMBS: usize, P> Copy for F2Ext<LIMBS, P> where P: BinaryIrreducible<LIMBS> {}

impl<const LIMBS: usize, P> PartialEq for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<const LIMBS: usize, P> Eq for F2Ext<LIMBS, P> where P: BinaryIrreducible<LIMBS> {}

impl<const LIMBS: usize, P> core::fmt::Debug for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "F2Ext({:?})", self.value)
    }
}

// ---------------------------------------------------------------------------
// Pretty Display — polynomial form over F_2
// ---------------------------------------------------------------------------
//
// Shows the element as a sum of powers of x, e.g.  `x^7 + x^3 + x + 1`.
// The zero element is printed as `0`.

impl<const LIMBS: usize, P> core::fmt::Display for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let m = P::degree();
        let mut terms = Vec::new();
        let words = self.value.to_words();

        // Collect set bits from high to low for descending-degree display.
        for i in (0..m).rev() {
            let word = words[i / 64];
            let bit = (word >> (i % 64)) & 1;
            if bit == 1 {
                match i {
                    0 => terms.push("1".to_string()),
                    1 => terms.push("x".to_string()),
                    _ => terms.push(format!("x^{i}")),
                }
            }
        }

        if terms.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", terms.join(" + "))
        }
    }
}

// ---------------------------------------------------------------------------
// CtOption functionalities
// ---------------------------------------------------------------------------

impl<const LIMBS: usize, P> Default for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn default() -> Self {
        Self::zero()
    }
}

impl<const LIMBS: usize, P> ConditionallySelectable for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self::new(Uint::<LIMBS>::conditional_select(
            &a.value, &b.value, choice,
        ))
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.value.conditional_assign(&other.value, choice)
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        Uint::<LIMBS>::conditional_swap(&mut a.value, &mut b.value, choice)
    }
}

impl<const LIMBS: usize, P> ConstantTimeEq for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        Uint::<LIMBS>::ct_eq(&self.value, &other.value)
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        Uint::<LIMBS>::ct_ne(&self.value, &other.value)
    }
}

// ---------------------------------------------------------------------------
// Private helper functions
// ---------------------------------------------------------------------------

fn reduce<const LIMBS: usize, P>(mut a: Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>,
{
    let modulus = P::modulus();
    let m = P::degree();
    let nbits = Uint::<LIMBS>::BITS as usize;

    debug_assert!(m > 0);
    debug_assert!(m < nbits);

    // Fixed iteration count => no leakage from the current degree of `a`.
    for i in (m..nbits).rev() {
        let bit = a.bit(i as u32); // Choice: whether x^i is present
        let shifted = modulus << (i - m);
        let reduced = a ^ shifted;

        // Constant-time:
        // - keep `a` if bit == 0
        // - replace by `reduced` if bit == 1
        a = Uint::<LIMBS>::conditional_select(&a, &reduced, bit.into());
    }
    a
}

fn add_helper<const LIMBS: usize>(a: &Uint<LIMBS>, b: &Uint<LIMBS>) -> Uint<LIMBS> {
    a ^ b
}

// This implements the add-and-shift with immediate modular reduction algorithm
fn mul_helper<const LIMBS: usize, P>(a: &Uint<LIMBS>, b: &Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>,
{
    let m = P::degree();
    let modulus = P::modulus();

    let mut res = Uint::<LIMBS>::ZERO;
    let mut cur = a.clone(); // current term appearing in the sum

    for i in 0..m {
        let bit = b.bit(i as u32);

        // if bit(b, i) = 1 then res += cur (sum here is just xor-ing)
        let res_xor = res ^ cur;
        res = Uint::<LIMBS>::conditional_select(&res, &res_xor, bit.into());

        // We now multiply cur by x and reduce modulo P

        let top = cur.bit((m - 1) as u32);
        // top represents the bit of cur at index (m-1), which is the highest coefficient of cur mod modulus

        let shifted = cur << 1;
        // if top = 0 then the x^(m-1) term in cur was 0,
        // so the term x^m does not appear in shifted and we don't need
        // to simplify shifted using the equality x^m = lower.
        let reduced = shifted ^ modulus;
        // if top = 1 then the term x^(m-1) appears in cur,
        // so x^m appears in shifted, and we need to simplify this using the equality x^m = lower.
        // It suffices to add modulus, because this adds lower and kills the leading term.
        cur = Uint::<LIMBS>::conditional_select(&shifted, &reduced, top.into());
    }

    res
}

fn square_helper<const LIMBS: usize, P>(a: &Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>,
{
    let m = P::degree();
    let modulus = P::modulus();

    let mut res = Uint::<LIMBS>::ZERO;
    let mut cur = Uint::<LIMBS>::ONE; // cur = x^(2i) mod P, initially i = 0

    for i in 0..m {
        let bit = a.bit(i as u32);

        // if bit(a, i) = 1 then res += x^(2i) mod P
        let res_xor = res ^ cur;
        res = Uint::<LIMBS>::conditional_select(&res, &res_xor, bit.into());

        // multiply cur by x^2 and reduce mod P (we do so in 2 steps)
        let top1 = cur.bit((m - 1) as u32);
        let shifted1 = cur << 1;
        let reduced1 = shifted1 ^ modulus;
        cur = Uint::<LIMBS>::conditional_select(&shifted1, &reduced1, top1.into());

        let top2 = cur.bit((m - 1) as u32);
        let shifted2 = cur << 1;
        let reduced2 = shifted2 ^ modulus;
        cur = Uint::<LIMBS>::conditional_select(&shifted2, &reduced2, top2.into());
    }

    res
}

fn pow_2k_helper<const LIMBS: usize, P>(a: &Uint<LIMBS>, k: &usize) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>,
{
    let mut res = a.clone();
    for _ in 0..*k {
        res = square_helper::<LIMBS, P>(&res);
    }
    res
}

fn itoh_tsujii<const LIMBS: usize, P>(a: &Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>,
{
    // the degree m is public so its bits are public too!
    let m = P::degree();
    if m == 1 {
        return a.clone();
    }

    let mut beta = a.clone();
    let mut r = 1usize;

    let top = (m - 1).ilog2(); // number of bits of m-1

    for i in (0..top).rev() {
        let beta_frob = pow_2k_helper::<LIMBS, P>(&beta, &r); // computing beta_r^(2r)
        beta = mul_helper::<LIMBS, P>(&beta, &beta_frob); // computing beta_(2r) = beta_r * beta_r^(2r)
        r <<= 1; // doubling r

        if (((m - 1) >> i) & 1) == 1 {
            beta = mul_helper::<LIMBS, P>(&square_helper::<LIMBS, P>(&beta), a);
            // the current bit being 1, we compute beta_(2r+1) = a * beta_(2r)^2
            r += 1;
        }
    }
    square_helper::<LIMBS, P>(&beta)
}

fn sqrt_helper<const LIMBS: usize, P>(a: &Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>,
{
    // In the binary field F_{2^m}, we have sqrt(a) = a^(2^(m-1)) for any a (even a = 0)

    let m = P::degree();
    pow_2k_helper::<LIMBS, P>(&a, &(m - 1))
}

// ===========================================================================
// Operator overloads (delegate to the FieldOps methods below)
// ===========================================================================

impl<'a, 'b, const LIMBS: usize, P> Add<&'b F2Ext<LIMBS, P>> for &'a F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = F2Ext<LIMBS, P>;

    fn add(self, rhs: &'b F2Ext<LIMBS, P>) -> Self::Output {
        <F2Ext<LIMBS, P> as FieldOps>::add(self, rhs)
    }
}

impl<'a, 'b, const LIMBS: usize, P> Sub<&'b F2Ext<LIMBS, P>> for &'a F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = F2Ext<LIMBS, P>;

    fn sub(self, rhs: &'b F2Ext<LIMBS, P>) -> Self::Output {
        <F2Ext<LIMBS, P> as FieldOps>::sub(self, rhs)
    }
}

impl<'a, 'b, const LIMBS: usize, P> Mul<&'b F2Ext<LIMBS, P>> for &'a F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = F2Ext<LIMBS, P>;

    fn mul(self, rhs: &'b F2Ext<LIMBS, P>) -> Self::Output {
        <F2Ext<LIMBS, P> as FieldOps>::mul(self, rhs)
    }
}

impl<'a, const LIMBS: usize, P> Neg for &'a F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = F2Ext<LIMBS, P>;

    fn neg(self) -> Self::Output {
        <F2Ext<LIMBS, P> as FieldOps>::negate(self)
    }
}

// ===========================================================================
// FieldOps implementation
// ===========================================================================

impl<const LIMBS: usize, P> FieldOps for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    fn zero() -> Self {
        Self::from_uint(Uint::<LIMBS>::ZERO)
    }

    fn one() -> Self {
        Self::from_uint(Uint::<LIMBS>::ONE)
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
        // we have x = -x for any element x of a binary field
        *self
    }

    fn add(&self, rhs: &Self) -> Self {
        Self::new(add_helper(&self.value, &rhs.value))
    }

    fn sub(&self, rhs: &Self) -> Self {
        Self::new(add_helper(&self.value, &rhs.value))
    }

    fn mul(&self, rhs: &Self) -> Self {
        Self::new(mul_helper::<LIMBS, P>(&self.value, &rhs.value))
    }

    fn square(&self) -> Self {
        Self::new(square_helper::<LIMBS, P>(&self.value))
    }

    fn double(&self) -> Self {
        // 2 = 0 in binary fields!
        Self::zero()
    }

    fn invert(&self) -> CtOption<Self> {
        let is_invertible = !self.is_zero();
        CtOption::new(
            Self::new(itoh_tsujii::<LIMBS, P>(&self.value)),
            is_invertible,
        )
    }

    fn frobenius(&self) -> Self {
        self.square()
    }

    fn trace(&self) -> Self {
        let deg = P::degree();
        let mut result = self.clone();
        let mut conj = self.frobenius();
        for _ in 1..deg {
            result = <Self as FieldOps>::add(&result, &conj);
            conj = conj.frobenius();
        }
        result
    }

    fn norm(&self) -> Self {
        let deg = P::degree();
        let mut result = self.clone();
        let mut conj = self.frobenius();
        for _ in 1..deg {
            result = <Self as FieldOps>::mul(&result, &conj);
            conj = conj.frobenius();
        }
        result
    }

    fn sqrt(&self) -> CtOption<Self> {
        // In binary fields, everything is a square, so the Choice object always wraps 1u8
        let sqrt = Self::new(sqrt_helper::<LIMBS, P>(&self.value));
        CtOption::new(sqrt, Choice::from(1u8))
    }

    fn legendre(&self) -> i8 {
        let is_zero = self.is_zero();
        i8::conditional_select(&0i8, &1i8, is_zero)
    }

    fn characteristic() -> Vec<u64> {
        vec![2u64]
    }

    fn degree() -> u32 {
        P::degree() as u32
    }
}

// ---------------------------------------------------------------------------
// Cryptographically secure random sampling
// ---------------------------------------------------------------------------

impl<const LIMBS: usize, P> FieldRandom for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    /// Sample a uniformly random element of F_{2^m} using a CSPRNG.
    ///
    /// Fills a `Uint<LIMBS>` with random bytes, masks to `m` bits,
    /// and wraps via `F2Ext::new` (which reduces mod the irreducible).
    fn random(rng: &mut (impl rand::CryptoRng + rand::Rng)) -> Self {
        let m = P::degree();
        let mut words = [0u64; LIMBS];
        for w in words.iter_mut() {
            *w = rng.next_u64();
        }

        // Mask to m bits so we stay in the valid range [0, 2^m).
        let full_limbs = m / 64;
        let leftover = m % 64;

        // Zero out limbs beyond the ones we need.
        for w in words
            .iter_mut()
            .skip(full_limbs + if leftover > 0 { 1 } else { 0 })
        {
            *w = 0;
        }
        // Mask the partial top limb.
        if leftover > 0 && full_limbs < LIMBS {
            words[full_limbs] &= (1u64 << leftover) - 1;
        }

        Self::new(Uint::from_words(words))
    }
}

impl<const LIMBS: usize, P> FieldFromRepr for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Repr = Uint<LIMBS>;

    fn from_repr(x: Self::Repr) -> Self {
        Self::from_uint(x)
    }
}
