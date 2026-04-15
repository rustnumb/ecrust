#![doc = include_str!("../../katex-header.html")]
// Generic finary fields F_{2^m} = F_2[x] / (f(x))

use crate::field_ops::FieldOps;
use core::ops::{Add, Mul, Neg, Sub};
use crypto_bigint::Uint;
use std::marker::PhantomData;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// ---------------------------------------------------------------------------
// IrreduciblePoly — the only thing callers need to implement for a new field
// ---------------------------------------------------------------------------

pub trait BinaryIrreducible<const LIMBS: usize>: 'static {
    // Full polynomial bitmask, including the leading term x^m
    fn modulus() -> Uint<LIMBS>;

    // Degree m of the irreducible polynomial
    fn degree() -> usize;
}

// ---------------------------------------------------------------------------
// F2Ext — element of F_{2^M}
// ---------------------------------------------------------------------------
pub struct F2Ext<const LIMBS: usize, P>
where
    P: BinaryIrreducible<LIMBS>,
{
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
    pub fn new(x: Uint<LIMBS>) -> Self {
        Self {
            value: reduce::<LIMBS, P>(x),
            _phantom: PhantomData,
        }
    }

    pub fn from_u64(x: u64) -> Self {
        Self::new(Uint::from(x))
    }

    pub fn from_uint(x: Uint<LIMBS>) -> Self {
        Self::new(x)
    }

    pub fn as_uint(&self) -> Uint<LIMBS> {
        self.value
    }

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

impl<const LIMBS: usize, P> Add for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FieldOps::add(&self, &rhs)
    }
}

impl<const LIMBS: usize, P> Sub for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FieldOps::sub(&self, &rhs)
    }
}

impl<const LIMBS: usize, P> Mul for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FieldOps::mul(&self, &rhs)
    }
}

impl<const LIMBS: usize, P> Neg for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    type Output = Self;
    fn neg(self) -> Self {
        FieldOps::negate(&self)
    }
}

// ===========================================================================
// FieldOps implementation
// ===========================================================================

impl<const LIMBS: usize, P> FieldOps for F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>,
{
    // to be done: legendre, sqrt
    fn zero() -> Self {
        Self::from_uint(Uint::<LIMBS>::ZERO)
    }

    fn one() -> Self {
        Self::from_uint(Uint::<LIMBS>::ONE)
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
