//! Generic extension field  Fp^M = Fp[x] / (f(x)).
//!
//! # Overview
//!
//! Given a prime p and a monic irreducible polynomial f of degree M over Fp,
//! the extension field Fp^M = Fp[x]/(f(x)) is the set of polynomials of
//! degree < M with coefficients in Fp, where arithmetic is done modulo f.
//!
//! This module provides a single generic type [`FpExt<MOD, LIMBS, M, P>`] that
//! covers *any* such extension.  The irreducible polynomial is supplied via the
//! zero-size marker trait [`IrreduciblePoly`].
//!
//! # Representation
//!
//! An element `a ∈ Fp^M` is stored as exactly M base-field coefficients:
//!
//! ```text
//! coeffs = [a_0, a_1, ..., a_{M-1}]
//!        ↔  a_0 + a_1x + ... + a_{M-1}x^{M-1}
//! ```
//!
//! # Operations and costs
//!
//! | Operation   | Algorithm                        | Base-field cost    |
//! |-------------|----------------------------------|--------------------|
//! | Add / Sub   | Coefficient-wise                 | M  adds            |
//! | Negate      | Coefficient-wise                 | M  negs            |
//! | Double      | Coefficient-wise                 | M  doubles         |
//! | Multiply    | Schoolbook + reduction mod f     | M^2 muls + M^2 adds  |
//! | Square      | Same as multiply (self * self)   | M^2 muls + M^2 adds  |
//! |
//!      | Polynomial extended GCD          | O(M^2)              |
//! | Frobenius   | self^p  via square-and-multiply  | O(M^2 log p)        |
//! | Norm        | Product of M Galois conjugates   | O(M^3 log p)        |
//! | Trace       | Sum of M Galois conjugates       | O(M^2 log p)        |

use core::ops::{Add, Mul, Neg, Sub};
use std::marker::PhantomData;

use crate::field_ops::FieldOps;
use crate::fp_element::FpElement;
use crypto_bigint::{modular::ConstPrimeMontyParams, Uint};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// ===========================================================================
// IrreduciblePoly — the only thing callers need to implement for a new field
// ===========================================================================

/// Supplies the monic irreducible polynomial
///   `f(x) = x^M + c_{M-1}x^{M-1} + ... + c_1x + c_0`
/// that defines the extension  `Fp^M = Fp[x] / (f(x))`.
///
/// # Convention
///
/// `modulus()` returns the **non-leading** coefficients in ascending degree order:
/// `[c_0, c_1, ..., c_{M-1}]`.
/// The leading coefficient `1` (coefficient of `x^M`) is implicit.
///
/// # Example: f(x) = x^2 + 1 over F_19
/// ```ignore
/// struct MyPoly;
/// impl IrreduciblePoly<Fp19Mod, 1, 2> for MyPoly {
///     fn modulus() -> [FpElement<Fp19Mod, 1>; 2] {
///         [FpElement::one(), FpElement::zero()]  // c_0=1, c_1=0 \to x^2+1
///     }
/// }
/// ```
///
/// # Example: f(x) = x^3 − 2 over F_19  (i.e. x^3 + 17 since −2 \equiv 17 mod 19)
/// ```ignore
/// struct MyCubicPoly;
/// impl IrreduciblePoly<Fp19Mod, 1, 3> for MyCubicPoly {
///     fn modulus() -> [FpElement<Fp19Mod, 1>; 3] {
///         // f(x) = x^3 + 0x^2 + 0x + 17  ->  [17, 0, 0]
///         [FpElement::from_u64(17), FpElement::zero(), FpElement::zero()]
///     }
/// }
/// ```
pub trait IrreduciblePoly<MOD, const LIMBS: usize, const M: usize>: 'static
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    /// Non-leading coefficients `[c_0, c_1, ..., c_{M-1}]` of the defining polynomial.
    fn modulus() -> [FpElement<MOD, LIMBS>; M];
}

/*
===========================================================================
TonelliShanksConstants - The only other thing users have to
implement, sorry
===========================================================================
*/

/// Supplies the group order of the multiplicative group
///
/// # Example: F_(19^3)
/// ```ignore
/// impl TonelliShanksConstants<3> for MyOrder {
///     // Still only need 1 limb for 19^3
///     fn order() -> Unit<1> {
///         const TSCONSTS: Uint<1> = Uint::<1>::from_u64(6858);
///         ORDER
///     }
///     fn half_order() -> Unit<1> {
///         const ORDER: Uint<1> = Uint::<1>::from_u64(3429);
///         ORDER
///     }
/// }
/// ```
pub trait TonelliShanksConstants<MOD, const LIMBS: usize, const M: usize, const N: usize>:
    'static
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    // Write the order of the multiplicative group as
    // (p^M - 1) = 2^S * T where T is odd
    // Multiplicative group order p^M - 1
    // p^M - 1
    const ORDER: Uint<N>;
    // (p^M - 1) / 2
    const HALF_ORDER: Uint<N>;
    // Constant S
    const S: u64;
    // Constant 2^(S - 1)
    const TWOSM1: Uint<N>;
    // Constant T
    const T: Uint<N>;
    // Projenator exponent of the TS algorithm this is (T - 1) / 2
    const PROJENATOR_EXP: Uint<N>;
    // Root of unity TODO: implement in a way in which this is a const
    fn root_of_unity() -> [FpElement<MOD, LIMBS>; M];
}

// ===========================================================================
// FpExt — element of Fp^M
// ===========================================================================

/// An element of the extension field  `Fp^M = Fp[x] / (f(x))`.
///
/// `P` is a zero-size marker type implementing [`IrreduciblePoly`].
/// `M` is the extension degree (number of base-field coefficients stored).
/// `N` is the number limbs needed to store p^M
pub struct FpExt<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    /// Coefficients in ascending degree: `coeffs[i]` = coefficient of `x^i`.
    pub coeffs: [FpElement<MOD, LIMBS>; M],
    _phantom: PhantomData<P>,
    _order: PhantomData<TSCONSTS>,
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS>
    FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    /// Construct from a coefficient array `[a_0, ..., a_{M-1}]`.
    pub fn new(coeffs: [FpElement<MOD, LIMBS>; M]) -> Self {
        Self {
            coeffs,
            _phantom: PhantomData,
            _order: PhantomData,
        }
    }

    /// Embed a base-field element as `a + 0x + ... + 0x^{M-1}`.
    pub fn from_base(a: FpElement<MOD, LIMBS>) -> Self {
        let mut coeffs = std::array::from_fn(|_| FpElement::zero());
        coeffs[0] = a;
        Self::new(coeffs)
    }

    /// Return the coefficient of `x^i`.
    pub fn coeff(&self, i: usize) -> &FpElement<MOD, LIMBS> {
        &self.coeffs[i]
    }
}

// ---------------------------------------------------------------------------
// Clone / Copy / PartialEq / Eq / Debug
// (manual impls so we don't over-constrain the bounds the way #[derive] would)
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Clone
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    fn clone(&self) -> Self {
        Self {
            coeffs: self.coeffs.clone(),
            _phantom: PhantomData,
            _order: PhantomData,
        }
    }
}

// FpElement is Copy (derives it), so [FpElement; M] is Copy too.
impl<MOD, const LIMBS: usize, const M: usize, P, const N: usize, TSCONSTS> Copy
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
    [FpElement<MOD, LIMBS>; M]: Copy,
{
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> PartialEq
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    fn eq(&self, other: &Self) -> bool {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .all(|(a, b)| a == b)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Eq
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> std::fmt::Debug
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
    FpElement<MOD, LIMBS>: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FpExt{:?}", self.coeffs.as_ref())
    }
}

// ---------------------------------------------------------------------------
// CtOption functionalities
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Default
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    fn default() -> Self {
        Self::zero()
    }
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> ConditionallySelectable
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let mut res_coeffs = [FpElement::zero(); M];
        for i in 0..M {
            res_coeffs[i] = FpElement::conditional_select(&a.coeffs[i], &b.coeffs[i], choice);
        }
        Self::new(res_coeffs)
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        for i in 0..M {
            self.coeffs[i].conditional_assign(&other.coeffs[i], choice);
        }
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        for i in 0..M {
            FpElement::conditional_swap(&mut a.coeffs[i], &mut b.coeffs[i], choice);
        }
    }
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> ConstantTimeEq
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    fn ct_eq(&self, other: &Self) -> Choice {
        let mut acc = Choice::from(1u8);
        for i in 0..M {
            acc &= self.coeffs[i].ct_eq(&other.coeffs[i]);
        }
        acc
    }

    fn ct_ne(&self, other: &Self) -> Choice {
        !self.ct_eq(other)
    }
}

// ===========================================================================
// Private polynomial helpers
//
// All intermediate polynomial arithmetic uses Vec<FpElement> so that we
// never need arithmetic on const generic bounds (e.g. [T; M+M-1]), which is
// not yet stable in Rust.  The conversion back to [T; M] happens only at the
// boundary (poly_reduce).
// ===========================================================================

type Poly<MOD, const LIMBS: usize> = Vec<FpElement<MOD, LIMBS>>;

/// Coefficient-wise addition.  Output length = max(|a|, |b|).
#[allow(dead_code)]
fn poly_add<MOD, const LIMBS: usize>(
    a: &[FpElement<MOD, LIMBS>],
    b: &[FpElement<MOD, LIMBS>],
) -> Poly<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    let len = a.len().max(b.len());
    (0..len)
        .map(|i| {
            let ai = a.get(i).cloned().unwrap_or_else(|| FpElement::zero());
            let bi = b.get(i).cloned().unwrap_or_else(|| FpElement::zero());
            FieldOps::add(&ai, &bi)
        })
        .collect()
}

/// Coefficient-wise subtraction.
fn poly_sub<MOD, const LIMBS: usize>(
    a: &[FpElement<MOD, LIMBS>],
    b: &[FpElement<MOD, LIMBS>],
) -> Poly<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    let len = a.len().max(b.len());
    (0..len)
        .map(|i| {
            let ai = a.get(i).cloned().unwrap_or_else(|| FpElement::zero());
            let bi = b.get(i).cloned().unwrap_or_else(|| FpElement::zero());
            FieldOps::sub(&ai, &bi)
        })
        .collect()
}

/// Multiply every coefficient by a scalar.
fn poly_scale<MOD, const LIMBS: usize>(
    a: &[FpElement<MOD, LIMBS>],
    s: &FpElement<MOD, LIMBS>,
) -> Poly<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    a.iter().map(|ai| FieldOps::mul(ai, s)).collect()
}

/// Schoolbook polynomial multiplication.  Output degree = deg(a) + deg(b).
fn poly_mul<MOD, const LIMBS: usize>(
    a: &[FpElement<MOD, LIMBS>],
    b: &[FpElement<MOD, LIMBS>],
) -> Poly<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut out = vec![FpElement::zero(); a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            let t = FieldOps::mul(ai, bj);
            out[i + j] = FieldOps::add(&out[i + j], &t);
        }
    }
    out
}

/// Remove trailing zero coefficients (i.e. normalise to the canonical form
/// where the leading coefficient is non-zero, or the empty Vec for the zero poly).
fn poly_normalize<MOD, const LIMBS: usize>(mut a: Poly<MOD, LIMBS>) -> Poly<MOD, LIMBS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    while a.last().map_or(false, |c| bool::from(c.is_zero())) {
        a.pop();
    }
    a
}

/// Polynomial long division.  Returns `(quotient, remainder)` such that
/// `a = quotient * b + remainder`  with  `deg(remainder) < deg(b)`.
///
/// Requires b ≠ 0 and its leading coefficient invertible (always true over Fp).
fn poly_divmod<MOD, const LIMBS: usize>(
    a: &[FpElement<MOD, LIMBS>],
    b: &[FpElement<MOD, LIMBS>],
) -> (Poly<MOD, LIMBS>, Poly<MOD, LIMBS>)
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    let b = poly_normalize(b.to_vec());
    assert!(!b.is_empty(), "poly_divmod: divisor is zero");

    let mut rem = poly_normalize(a.to_vec());

    if rem.len() < b.len() {
        return (vec![], rem);
    }

    let b_lead_inv = b
        .last()
        .unwrap()
        .invert()
        .expect("poly_divmod: leading coefficient of b is not invertible");

    let out_len = rem.len() - b.len() + 1;
    let mut quot = vec![FpElement::zero(); out_len];

    // Process from the highest-degree term of the quotient down to degree 0.
    for i in (0..out_len).rev() {
        // The current highest remaining degree is i + deg(b).
        let rem_deg = i + b.len() - 1;
        let coeff = FieldOps::mul(&rem[rem_deg], &b_lead_inv);
        if bool::from(coeff.is_zero()) {
            continue;
        }
        quot[i] = coeff.clone();
        for (j, bj) in b.iter().enumerate() {
            let t = FieldOps::mul(&coeff, bj);
            rem[i + j] = FieldOps::sub(&rem[i + j], &t);
        }
    }

    (quot, poly_normalize(rem))
}

/// Reduce polynomial `a` modulo the irreducible `f(x) = x^M + Σ modulus[j]x^j`.
///
/// Uses the substitution rule:
///   `x^M \equiv −(modulus[0] + modulus[1]x + ... + modulus[M-1]x^{M-1})`
///
/// Sweeps from the highest degree of `a` down to `M`, eliminating one term
/// per step.  Each step costs M multiplications and M additions in Fp.
fn poly_reduce<MOD, const LIMBS: usize, const M: usize>(
    a: Vec<FpElement<MOD, LIMBS>>,
    modulus: &[FpElement<MOD, LIMBS>; M],
) -> [FpElement<MOD, LIMBS>; M]
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    let mut a = a;
    // Ensure we have at least M slots.
    while a.len() < M {
        a.push(FpElement::zero());
    }

    // Eliminate all terms of degree >= M, from the top down.
    for i in (M..a.len()).rev() {
        let lead = a[i].clone();
        if bool::from(lead.is_zero()) {
            continue;
        }
        // x^i = x^{i−M}  x^M  \equiv  −x^{i−M}  Σ_j modulus[j]x^j
        // → subtract leadmodulus[j] from position (i−M+j) for each j.
        for j in 0..M {
            let t = FieldOps::mul(&lead, &modulus[j]);
            a[i - M + j] = FieldOps::sub(&a[i - M + j], &t);
        }
        a[i] = FpElement::zero();
    }

    // Copy the first M coefficients into a fixed-size array.
    std::array::from_fn(|i| {
        if i < a.len() {
            a[i].clone()
        } else {
            FpElement::zero()
        }
    })
}

/// Polynomial extended GCD over Fp.
///
/// Returns `(g, s)` satisfying  `as \equiv g  (mod b)`,  where `g = gcd(a, b)`.
///
/// When `b` is irreducible and `a` is not a multiple of `b`, `g` is a nonzero
/// constant (units are the GCD of coprime polynomials over a field), so
/// `a⁻¹ \equiv s  g⁻¹  (mod b)`.
///
/// Uses the standard iterative Euclidean algorithm:
/// ```text
/// (r₀, s₀) = (a, 1)
/// (r₁, s₁) = (b, 0)
/// while r₁ ≠ 0:
///     q         = r₀ / r₁
///     (r₀, r₁) = (r₁, r₀ − qr₁)
///     (s₀, s₁) = (s₁, s₀ − qs₁)
/// return (r₀, s₀)
/// ```
fn poly_ext_gcd<MOD, const LIMBS: usize>(
    a: Poly<MOD, LIMBS>,
    b: Poly<MOD, LIMBS>,
) -> (Poly<MOD, LIMBS>, Poly<MOD, LIMBS>)
// (gcd, s)
where
    MOD: ConstPrimeMontyParams<LIMBS>,
{
    let mut r0 = poly_normalize(a);
    let mut r1 = poly_normalize(b);
    let mut s0: Poly<MOD, LIMBS> = vec![FpElement::one()]; // 1
    let mut s1: Poly<MOD, LIMBS> = vec![]; // 0

    while !r1.is_empty() {
        let (q, r) = poly_divmod(&r0, &r1);
        // s_new = s0 − qs1
        let q_s1 = poly_mul(&q, &s1);
        let s_new = poly_normalize(poly_sub(&s0, &q_s1));
        s0 = s1;
        s1 = s_new;
        r0 = r1;
        r1 = poly_normalize(r);
    }

    (r0, s0)
}

// ===========================================================================
// Operator overloads (delegate to the FieldOps methods below)
// ===========================================================================

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Add
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FieldOps::add(&self, &rhs)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Sub
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FieldOps::sub(&self, &rhs)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Mul
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FieldOps::mul(&self, &rhs)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> Neg
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    type Output = Self;
    fn neg(self) -> Self {
        FieldOps::negate(&self)
    }
}

// ===========================================================================
// Helper functions for Tonelli-Shanks algorithm
// ===========================================================================

/// The loop part of the Tonelli--Shanks algorithm
///
/// Takes as input a (reference to an) element `x` of FpM and the
/// projenator `w = x^((t-1)/2)` for `x`. The function then modifies
/// `x` and return the final value sqrt(x) if `x` is a quadratic
/// residue in FpM
///
/// # Arguments
///
/// * `x` - The value to take the sqrt of (type: &mut FpExt)
/// * `w` - The projenator for x (type: &FpExt)
///
/// # Note
///
/// Most of this is directly taken from the Fp case at
/// https://github.com/SamFrengley/ff-sqrtratio/blob/8e5aa6f934f32d9b9cff56177d9943a2effcd390/ff_derive/src/lib.rs
/// under an MIT liscence.
fn ts_loop<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS>(
    x: &mut FpExt<MOD, LIMBS, M, N, P, TSCONSTS>,
    w: &FpExt<MOD, LIMBS, M, N, P, TSCONSTS>,
) where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    let mut v = TSCONSTS::S as u32;
    *x = x.mul(*w);
    let mut b = x.mul(*w);
    let mut z = FpExt::new(TSCONSTS::root_of_unity());

    for max_v in (1..=TSCONSTS::S as u32).rev() {
        let mut k = 1u32;
        let mut tmp = b.square();
        let mut j_less_than_v: Choice = 1.into();

        for j in 2..max_v {
            let tmp_is_one = tmp.ct_eq(&FpExt::one());
            let squared = FpExt::conditional_select(&tmp, &z, tmp_is_one).square();
            tmp = FpExt::conditional_select(&squared, &tmp, tmp_is_one);
            let new_z = FpExt::conditional_select(&z, &squared, tmp_is_one);

            j_less_than_v &= !j.ct_eq(&v);

            k = u32::conditional_select(&j, &k, tmp_is_one);
            z = FpExt::conditional_select(&z, &new_z, j_less_than_v);
        }

        let result = x.mul(z);
        *x = FpExt::conditional_select(&result, x, b.ct_eq(&FpExt::one()));
        z = z.square();
        b = b.mul(z);
        v = k;
    }
}

// ===========================================================================
// FieldOps implementation
// ===========================================================================

impl<MOD, const LIMBS: usize, const M: usize, const N: usize, P, TSCONSTS> FieldOps
    for FpExt<MOD, LIMBS, M, N, P, TSCONSTS>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    TSCONSTS: TonelliShanksConstants<MOD, LIMBS, M, N>,
{
    // --- Identity elements --------------------------------------------------

    fn zero() -> Self {
        Self::new(std::array::from_fn(|_| FpElement::zero()))
    }

    fn one() -> Self {
        let mut c: [FpElement<MOD, LIMBS>; M] = std::array::from_fn(|_| FpElement::zero());
        c[0] = FpElement::one();
        Self::new(c)
    }

    // --- Predicates ---------------------------------------------------------

    fn is_zero(&self) -> Choice {
        self.ct_eq(&Self::zero())
    }

    fn is_one(&self) -> Choice {
        Self::ct_eq(self, &Self::one())
    }

    // --- Core arithmetic ----------------------------------------------------

    fn negate(&self) -> Self {
        Self::new(std::array::from_fn(|i| self.coeffs[i].negate()))
    }

    fn add(&self, rhs: &Self) -> Self {
        Self::new(std::array::from_fn(|i| {
            FieldOps::add(&self.coeffs[i], &rhs.coeffs[i])
        }))
    }

    fn sub(&self, rhs: &Self) -> Self {
        Self::new(std::array::from_fn(|i| {
            FieldOps::sub(&self.coeffs[i], &rhs.coeffs[i])
        }))
    }

    /// Schoolbook multiplication followed by reduction modulo f(x).
    ///
    /// Product has degree ≤ 2M−2, then each high-degree term is replaced using
    /// `x^M \equiv −Σ modulus[j]x^j` until all degrees are below M.
    fn mul(&self, rhs: &Self) -> Self {
        let product = poly_mul(&self.coeffs, &rhs.coeffs);
        Self::new(poly_reduce(product, &P::modulus()))
    }

    fn square(&self) -> Self {
        // Calls our own FieldOps::mul, not Mul::mul (operator), to avoid the
        // supertrait ambiguity documented in field_ops.rs.
        <Self as FieldOps>::mul(self, self)
    }

    fn double(&self) -> Self {
        Self::new(std::array::from_fn(|i| self.coeffs[i].double()))
    }

    // --- Inversion ----------------------------------------------------------

    /// Inversion via polynomial extended GCD.
    ///
    /// Finds `s(x)` such that `self(x)s(x) \equiv 1  (mod f(x))` by computing
    /// `gcd(self, f) = g` (a nonzero constant if self ≠ 0) and setting
    /// `self⁻¹ = sg⁻¹ mod f`.
    fn invert(&self) -> CtOption<Self> {
        let is_invertible = !self.is_zero();

        // Build the full irreducible polynomial as a Vec:
        // f = [c_0, c_1, ..., c_{M-1}, 1]  (coefficients in ascending degree)
        let mut f: Poly<MOD, LIMBS> = P::modulus().iter().cloned().collect();
        f.push(FpElement::one()); // leading x^M term

        let a: Poly<MOD, LIMBS> = self.coeffs.iter().cloned().collect();
        let (gcd, s) = poly_ext_gcd(a, f);

        // gcd is a nonzero constant in Fp (since f is irreducible and a ≠ 0).
        let g0 = gcd.into_iter().next().unwrap_or(FpElement::zero()); // the constant term = gcd value
        let g_inv = g0.invert().unwrap_or(FpElement::zero()); // invert in Fp

        // self⁻¹ = s(x)  g⁻¹  reduced mod f
        let s_scaled = poly_scale(&s, &g_inv);
        CtOption::new(
            Self::new(poly_reduce(s_scaled, &P::modulus())),
            is_invertible,
        )
    }

    // --- Frobenius ----------------------------------------------------------

    /// `φ_p(a) = a^p` — the p-power Frobenius endomorphism.
    ///
    /// Computed via square-and-multiply using the characteristic p retrieved
    /// from the base field.
    fn frobenius(&self) -> Self {
        let p = FpElement::<MOD, LIMBS>::characteristic();
        <Self as FieldOps>::pow(self, &p)
    }

    // --- Norm and trace -----------------------------------------------------

    /// `N_{Fp^M/Fp}(a) = ∏_{i=0}^{M-1} φ_p^i(a)` — product of all Galois conjugates.
    ///
    /// The result lies in Fp (all higher coefficients are 0), but is returned
    /// embedded in Fp^M for uniformity with the [`FieldOps`] signature.
    fn norm(&self) -> Self {
        let mut result = self.clone();
        let mut conj = self.frobenius();
        for _ in 1..M {
            result = <Self as FieldOps>::mul(&result, &conj);
            conj = conj.frobenius();
        }
        result
    }

    /// `Tr_{Fp^M/Fp}(a) = Σ_{i=0}^{M-1} φ_p^i(a)` — sum of all Galois conjugates.
    ///
    /// Like `norm`, the result lies in Fp but is returned embedded in Fp^M.
    fn trace(&self) -> Self {
        let mut result = self.clone();
        let mut conj = self.frobenius();
        for _ in 1..M {
            result = <Self as FieldOps>::add(&result, &conj);
            conj = conj.frobenius();
        }
        result
    }

    // --- Square root --------------------------------------------------------

    /// Tonelli--Shanks squareroot algorithm
    ///
    /// Implementation of the Tonelli--Shanks square root algorithm. Requires
    /// only a factorisation as $p^M - 1 = 2^K * N$ so can compute this at
    /// compile time by truncating zeros.
    ///
    /// # Arguments
    ///
    /// * `&self` - An element of Fp^M (type: &self)
    ///
    /// # Returns
    ///
    /// `self^(1/2)` a choice of squareroot (type: Self)
    fn sqrt(&self) -> CtOption<Self> {
        let mut x = *self;
        let exp = TSCONSTS::PROJENATOR_EXP;
        let exp_limbs = exp.as_limbs().map(|limb| limb.0);
        // TODO: generate the addition chain for this specific constant
        // this is constant time since exp_limbs is always(!) the same
        let w = x.pow_vartime(&exp_limbs);
        ts_loop(&mut x, &w);
        CtOption::new(
            x,
            x.mul(x).ct_eq(self), // Only return Some if it's the square root.
        )
    }

    /// Inverse and sqrt in one exponentiation
    ///
    /// Computes the inverse and squareroot of `self` in one
    /// exponentiation using the tricks in Scott's article
    ///
    /// # Arguments
    ///
    /// * `&self` - An element of Fp^M (type: &self)
    ///
    /// # Returns
    ///
    /// `(myinv, mysqrt)` which is `self.invert()` and `self.sqrt()`
    fn inverse_and_sqrt(&self) -> (CtOption<Self>, CtOption<Self>) {
        let is_invertible = !self.is_zero();

        let mut mysqrt = *self;
        let exp = TSCONSTS::PROJENATOR_EXP;
        let exp_limbs = exp.as_limbs().map(|limb| limb.0);
        // TODO: generate the addition chain for this specific constant
        // this is constant time since exp_limbs is always(!) the same
        let w = mysqrt.pow_vartime(&exp_limbs);
        ts_loop(&mut mysqrt, &w);

        let e0 = TSCONSTS::TWOSM1;
        let e0_limbs = e0.as_limbs().map(|limb| limb.0);
        let e1 = e0.sub(Uint::<N>::from_u64(1));
        let e1_limbs = e1.as_limbs().map(|limb| limb.0);

        // Compute x^(2^(S - 1) - 1) * (x * w^4)^(2^(S - 1))
        let myinv = self
            .pow_vartime(&e1_limbs)
            .mul((self.mul(&w.pow_vartime(&[4 as u64]))).pow_vartime(&e0_limbs));

        (
            CtOption::new(myinv, is_invertible),
            CtOption::new(
                mysqrt,
                mysqrt.mul(mysqrt).ct_eq(self), // Only return Some if it's the square root.
            ),
        )
    }

    /// Inverse of squareroot of `self` in 1 exponentiation
    ///
    /// Computes 1/sqrt(self) using the trick from Mike Scott's
    /// "Tricks of the trade" article Section 2
    /// https://eprint.iacr.org/2020/1497
    ///
    /// # Arguments
    ///
    /// * `&self` - Description of &self (type: self)
    ///
    /// # Returns
    ///
    /// The inverse of the squareroot of `self` (type: CtOption<Self>)
    fn inv_sqrt(&self) -> CtOption<Self> {
        let (inv, sqrt) = self.inverse_and_sqrt();
        inv.and_then(|a| sqrt.map(|b| a * b))
    }

    /// Inverse of `self` and squareroot of `rhs` in 1 exponentiation
    ///
    /// Computes `1/self` and `rhs.sqrt()` simulaineously using the
    /// trick from Mike Scott's "Tricks of the trade" article Section
    /// 2 https://eprint.iacr.org/2020/1497
    ///
    /// # Returns
    ///
    /// The inverse of `self` and square root fo `rhs`. Theq former is
    /// none if and only if `self` is nonzero and the latter is not
    /// none if and only if there exists a squareroot of `rhs` in FpM
    /// (type: (CtOption<Self>, CtOption<self>))
    fn invertme_sqrtother(&self, rhs: &Self) -> (CtOption<Self>, CtOption<Self>) {
        let is_invertible = !self.is_zero();

        let x = self.mul(self).mul(*rhs);
        let (myinv, mysqrt) = x.inverse_and_sqrt();

        let myinv_value = myinv.unwrap_or(Self::zero());
        let mysqrt_value = mysqrt.unwrap_or(Self::zero());
        let inv_value = self.mul(rhs).mul(myinv_value);
        let sqrt_value = inv_value.mul(mysqrt_value);

        let inv_is_some = is_invertible & myinv.is_some();

        let sqrt_is_some = inv_is_some & mysqrt.is_some() & (sqrt_value.mul(sqrt_value)).ct_eq(rhs);

        (
            CtOption::new(inv_value, inv_is_some),
            CtOption::new(sqrt_value, sqrt_is_some),
        )
    }

    /// Computes the squareroot of a ratio `self/rhs`
    ///
    /// Computes `sqrt(self/rhs)` in one exponentiation using the
    /// trick from Mike Scott's "Tricks of the trade" article Section
    /// 2 https://eprint.iacr.org/2020/1497
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
        let x = self.mul(&self.mul(self)).mul(*rhs);
        let (myinv, mysqrt) = x.inverse_and_sqrt();

        let myinv_value = myinv.unwrap_or(Self::zero());
        let mysqrt_value = mysqrt.unwrap_or(Self::zero());
        let ans_value = self.mul(self).mul(myinv_value).mul(mysqrt_value);

        let inv_is_some = myinv.is_some();
        let ans_is_some =
            inv_is_some & mysqrt.is_some() & (mysqrt_value.mul(mysqrt_value)).ct_eq(&x);

        CtOption::new(ans_value, ans_is_some)
    }

    /// a is a QR in Fp^M iff a^{(p^M-1)/2} = 1.
    ///
    /// Implements the "Legendre symbol" which is 1 if and only if we
    /// have a quadratic residue in FpM
    /// WARNING: Not constant time if `self` is zero
    ///
    /// # Arguments
    ///
    /// * `&self` - Element of FpM (type: self)
    ///
    /// # Returns
    ///
    /// Either `0` if `&self` is `0`, `1` if `&self` is a QR or `-1` if
    /// `&self` is not a QR. (type: i8)
    fn legendre(&self) -> i8 {
        let exp = TSCONSTS::HALF_ORDER;
        let exp_limbs = exp.as_limbs().map(|limb| limb.0);
        let symb = self.pow(&exp_limbs); // note, this is constant time since exp is constant

        let ret = i8::conditional_select(&-1, &1, symb.is_one());
        i8::conditional_select(&0, &ret, !symb.is_zero())
    }

    // --- Utilities ----------------------------------------------------------

    fn characteristic() -> Vec<u64> {
        // Same prime p as the base field.
        FpElement::<MOD, LIMBS>::characteristic()
    }

    /// Degree of Fp^M over the prime subfield Fp.
    fn degree() -> u32 {
        M as u32
    }
}
