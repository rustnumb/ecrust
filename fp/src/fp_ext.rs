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
//! | Invert      | Polynomial extended GCD          | O(M^2)              |
//! | Frobenius   | self^p  via square-and-multiply  | O(M^2 log p)        |
//! | Norm        | Product of M Galois conjugates   | O(M^3 log p)        |
//! | Trace       | Sum of M Galois conjugates       | O(M^2 log p)        |

use core::ops::{Add, Mul, Neg, Sub};
use std::marker::PhantomData;

use crypto_bigint::modular::ConstPrimeMontyParams;

use crate::field_ops::FieldOps;
use crate::fp_element::FpElement;

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

// ===========================================================================
// FpExt — element of Fp^M
// ===========================================================================

/// An element of the extension field  `Fp^M = Fp[x] / (f(x))`.
///
/// `P` is a zero-size marker type implementing [`IrreduciblePoly`].
/// `M` is the extension degree (number of base-field coefficients stored).
pub struct FpExt<MOD, const LIMBS: usize, const M: usize, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    /// Coefficients in ascending degree: `coeffs[i]` = coefficient of `x^i`.
    pub coeffs: [FpElement<MOD, LIMBS>; M],
    _phantom: PhantomData<P>,
}

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<MOD, const LIMBS: usize, const M: usize, P> FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    /// Construct from a coefficient array `[a_0, ..., a_{M-1}]`.
    pub fn new(coeffs: [FpElement<MOD, LIMBS>; M]) -> Self {
        Self {
            coeffs,
            _phantom: PhantomData,
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

impl<MOD, const LIMBS: usize, const M: usize, P> Clone for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    fn clone(&self) -> Self {
        Self {
            coeffs: self.coeffs.clone(),
            _phantom: PhantomData,
        }
    }
}

// FpElement is Copy (derives it), so [FpElement; M] is Copy too.
impl<MOD, const LIMBS: usize, const M: usize, P> Copy for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    [FpElement<MOD, LIMBS>; M]: Copy,
{
}

impl<MOD, const LIMBS: usize, const M: usize, P> PartialEq for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    fn eq(&self, other: &Self) -> bool {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .all(|(a, b)| a == b)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, P> Eq for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
}

impl<MOD, const LIMBS: usize, const M: usize, P> std::fmt::Debug for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
    FpElement<MOD, LIMBS>: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FpExt{:?}", self.coeffs.as_ref())
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
    while a.last().map_or(false, |c| c.is_zero()) {
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
        if coeff.is_zero() {
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
        if lead.is_zero() {
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

impl<MOD, const LIMBS: usize, const M: usize, P> Add for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FieldOps::add(&self, &rhs)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, P> Sub for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FieldOps::sub(&self, &rhs)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, P> Mul for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FieldOps::mul(&self, &rhs)
    }
}

impl<MOD, const LIMBS: usize, const M: usize, P> Neg for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
{
    type Output = Self;
    fn neg(self) -> Self {
        FieldOps::negate(&self)
    }
}

// ===========================================================================
// FieldOps implementation
// ===========================================================================

impl<MOD, const LIMBS: usize, const M: usize, P> FieldOps for FpExt<MOD, LIMBS, M, P>
where
    MOD: ConstPrimeMontyParams<LIMBS>,
    P: IrreduciblePoly<MOD, LIMBS, M>,
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

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|c| c.is_zero())
    }

    fn is_one(&self) -> bool {
        self.coeffs[0].is_one() && self.coeffs[1..].iter().all(|c| c.is_zero())
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
    fn invert(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Build the full irreducible polynomial as a Vec:
        // f = [c_0, c_1, ..., c_{M-1}, 1]  (coefficients in ascending degree)
        let mut f: Poly<MOD, LIMBS> = P::modulus().iter().cloned().collect();
        f.push(FpElement::one()); // leading x^M term

        let a: Poly<MOD, LIMBS> = self.coeffs.iter().cloned().collect();
        let (gcd, s) = poly_ext_gcd(a, f);

        // gcd is a nonzero constant in Fp (since f is irreducible and a ≠ 0).
        let g0 = gcd.into_iter().next()?; // the constant term = gcd value
        let g_inv = g0.invert()?; // invert in Fp

        // self⁻¹ = s(x)  g⁻¹  reduced mod f
        let s_scaled = poly_scale(&s, &g_inv);
        Some(Self::new(poly_reduce(s_scaled, &P::modulus())))
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
    fn sqrt(&self) -> Self {
        let p = FpElement::<MOD, LIMBS>::characteristic();
        let q = (&p, M);
        assert!();
    }

    fn legendre(&self) -> i8 {
        // a is a QR in Fp^M iff a^{(p^M-1)/2} = 1.
        // Placeholder: implement per-field when needed.
        todo!("FpExt::legendre — implement for your specific (p, M)")
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
