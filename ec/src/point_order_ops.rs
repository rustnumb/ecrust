//! Order computations for points on elliptic curves.
//!
//! # Problem
//!
//! Given a curve $E/\mathbb{F}_q$ and a point $P \in E(\mathbb{F}_q)$, we want
//! to compute the **order** of $P$: the smallest integer $n \ge 1$ such that
//!
//! $$
//! [n] P = O.
//! $$
//!
//! Following Sutherland (2007), we expose two regimes:
//!
//! | Regime                | Input given                                                                          | Cost                                                             |
//! |-----------------------|--------------------------------------------------------------------------------------|------------------------------------------------------------------|
//! | Known group order     | $N = \#E(\mathbb{F}_q)$ and its factorization                                        | $O(k \log N)$ scalar muls, where $k = \omega(N)$                 |
//! | Generic (Hasse bound) | Interval $[\ell_{lo}, \ell_{hi}]$ with $\lvert P \rvert \in [\ell_{lo}, \ell_{hi}]$  | $\widetilde O\!\left(\sqrt{\ell_{hi} - \ell_{lo}}\right)$        |
//!
//! # When to use which
//!
//! - If you already know $N = \#E(\mathbb{F}_q)$ (e.g. via SEA in the isogeny
//!   crate), use [`order_from_group_order`]. This is the usual cryptographic
//!   case and runs in time polynomial in $\log q$.
//! - Otherwise, bracket $\lvert P \rvert$ by the Hasse interval
//!   $[q + 1 - 2\sqrt q,\; q + 1 + 2\sqrt q]$ and call
//!   [`order_in_interval`].
//!
//! # On the complexity claim
//!
//! The present [`order_in_interval`] gives the textbook $\widetilde O(\sqrt W)$
//! where $W = \ell_{hi} - \ell_{lo}$ is the Hasse width (so $W \in O(\sqrt q)$
//! and the total cost is $\widetilde O(q^{1/4})$).
//!
//! Sutherland (2007, §5) shows this can be refined to
//! $\widetilde O\!\left(\sqrt{\ell_{\max}(\lvert P \rvert)}\right)$ by
//! peeling off the smooth part of $\lvert P \rvert$ before running BSGS,
//! which is strictly faster whenever $\lvert P \rvert$ has any nontrivial
//! smooth factor. That optimisation is non-trivial to get right — a naive
//! "test $\ell$-divisibility on the residual" does not work from a Hasse
//! interval alone — and is deferred to a follow-up implementation. See
//! the `TODO(sutherland-phase-1)` note at the bottom of this module.
//!
//! # Model coverage
//!
//! The algorithms are written against [`PointOps`] and [`PointAdd`] so they
//! work for every curve model in this crate that supports full group
//! addition: Weierstrass, Edwards, Hessian (and twisted Hessian), Jacobi
//! quartic, Jacobi intersection, and Legendre.
//!
//! [`order_in_interval`] additionally requires `P: Eq + Hash` for the
//! baby-step / giant-step table.  Deriving `Hash` on the point structs is a
//! one-line change wherever the underlying field element type is `Hash`
//! (which is the case for `Uint`-backed `FpElement`).
//!
//! Montgomery / Kummer points (x-only) do **not** implement [`PointAdd`] and
//! need a dedicated x-only variant using the pseudo-addition ladder; that
//! specialisation is left for a follow-up.

use crate::num_utils::{product_of_factors, trial_factor};
use crate::point_ops::{PointAdd, PointOps};
use crypto_bigint::{NonZero, Uint};
use std::collections::HashMap;
use std::hash::Hash;

// ---------------------------------------------------------------------------
// Case A:  known group order    (fast path,  O(log q))
// ---------------------------------------------------------------------------

/// Given a point `P`, the group order `N = #E(𝔽_q)`, and the prime
/// factorization of `N` as `[(p₁, e₁), …, (p_k, e_k)]`, return the exact
/// order of `P`.
///
/// # Algorithm
///
/// Since $\lvert P \rvert \mid N$ (Lagrange), start from $n \leftarrow N$
/// and, for each prime factor $p_i$ of $N$, strip as many copies of $p_i$
/// from $n$ as possible while still satisfying $[n] P = O$:
///
/// ```text
/// n ← N
/// for each (pᵢ, eᵢ) in factors:
///     for j = 1 .. eᵢ:
///         if [n / pᵢ] P = O:
///             n ← n / pᵢ
///         else:
///             break
/// return n
/// ```
///
/// The final `n` is the exact order of `P`.
///
/// # Cost
///
/// $O(k  log N)$ scalar multiplications, where $k = \omega(N)$ is the number
/// of distinct prime factors of `N`.
///
///
pub fn order_from_group_order<P, const L: usize>(
    point: &P,
    curve: &P::Curve,
    n_group: &Uint<L>,
    factors: &[(Uint<L>, u32)],
) -> Uint<L>
where
    P: PointOps,
{
    debug_assert_eq!(product_of_factors::<L>(factors), *n_group);

    let mut n = *n_group;

    //  For each (pᵢ, eᵢ), peel copies of pᵢ off n while [n/pᵢ] P = O.
    //  At most eᵢ successful peels per prime  ⇒ O(∑ eᵢ) = O(log N) scalar
    //  multiplications total.
    for (p, _e) in factors {
        let p_nz = NonZero::new(*p).expect("prime factor must be nonzero");
        loop {
            //  q = n / p.  If p does not divide n we stop (shouldn't
            //  happen for a prime factor of N, but guards bad input).
            let (q, rem) = n.div_rem(&p_nz);
            if !bool::from(rem.is_zero()) {
                break;
            }
            let q_point = point.scalar_mul(q.as_words(), curve);
            if q_point.is_identity() {
                n = q;
            } else {
                break;
            }
        }
    }

    n
}

// ---------------------------------------------------------------------------
// Case B:  generic order in a bracket   (BSGS, O(√W))
// ---------------------------------------------------------------------------

/// Compute the order of `P` knowing only a bracket $[\ell_{lo}, \ell_{hi}]$
/// that contains it.
///
/// For an elliptic curve over $\mathbb{F}_q$ the natural bracket is the
/// **Hasse interval**:
///
/// $$
/// \lvert P \rvert \in [\,q + 1 - 2\sqrt q,\; q + 1 + 2\sqrt q\,].
/// $$
///
/// # Algorithm (baby-step / giant-step)
///
/// Let $W = \ell_{hi} - \ell_{lo}$ and $w = \lceil \sqrt W \rceil$.
/// We look for a pair $(i, j)$ with $0 \le i, j \le w$ satisfying
///
/// $$
/// [\ell_{lo} + i \cdot w - j] P = O,
/// $$
///
/// which is equivalent to the collision
///
/// $$
/// [\ell_{lo} + i \cdot w]\, P \;=\; [j]\, P
/// $$
///
/// in the group. The baby table $\{[j] P : 0 \le j \le w\}$ is built
/// once; giant steps walk from $[\ell_{lo}] P$ in strides of $[w] P$
/// until a collision is found.
///
/// ```text
/// w ← ⌈√(\ell_hi − \ell_lo)⌉
/// table ← { [j] P ↦ j : 0 ≤ j ≤ w }
/// Q ← [\ell_lo] P
/// for i = 0 .. w:
///     if Q ∈ table with value j:
///         return  reduce_order( \ell_lo + i·w − j )
///     Q ← Q + [w] P
/// ```
///
/// The returned value $c = \ell_{lo} + i \cdot w - j$ need not itself be
/// the order; it is only guaranteed to be a multiple of it inside the
/// bracket. We finish by factoring $c$ and running the Case-A peel.
///
/// # Cost
///
/// `O(\sqrtW)` point additions and a hash table of the same size for the
/// BSGS itself. Phase 3 additionally trial-factors the resulting
/// witness `c ≤ hi`, which costs `O(\sqrtc · M)` big-int operations —
/// practical only when `hi` fits in roughly 80 bits.
///
/// In cryptographic settings ($q \ge 2^{128}$) you should **always**
/// prefer [`order_from_group_order`], feeding in a factorization you
/// obtained separately (e.g. from SEA). `order_in_interval` is
/// intended for testing, for small-field experimentation, and as the
/// raw building block that the future Sutherland smoothness
/// optimisation will wrap.
///
/// # Parameters
///
/// - `point` — the point `P` whose order we want
/// - `curve` — the curve it lives on
/// - `lo, hi` — bracket containing `|P|`
///
/// # Returns
///
/// The exact order `|P|`.
///
/// # Panics
///
/// - If `[k] P ≠ O` for every `k ∈ [lo, hi]`, i.e. the bracket is wrong.
/// - If the BSGS table size would exceed `u64::MAX` (unreachable in
///   practice: for any Hasse interval of a curve over a field fitting in
///   `Uint<L>`, $\sqrt W$ fits comfortably in 64 bits).
pub fn order_in_interval<P, const L: usize>(
    point: &P,
    curve: &P::Curve,
    lo: &Uint<L>,
    hi: &Uint<L>,
) -> Uint<L>
where
    P: PointAdd + Eq + Hash,
{
    //  Giant stride  w = ⌈\sqrt(hi − lo)⌉ .
    let width = hi.wrapping_sub(lo);
    let w_floor = width.floor_sqrt_vartime();
    let w_uint = if w_floor.wrapping_mul(&w_floor) == width {
        w_floor
    } else {
        w_floor.wrapping_add(&Uint::<L>::ONE)
    };

    //  BSGS table size fits in u64 for any realistic q — width ≲ 4\sqrtq
    //  so \sqrtwidth ≲ 2 · q^{1/4}, which is < 2^64 for q < 2^254.
    debug_assert!(
        w_uint.as_words()[1..].iter().all(|&w| w == 0),
        "BSGS table size would exceed u64; bracket is unrealistically large"
    );
    let w: u64 = w_uint.as_words()[0];

    //  BSGS: find candidate multiple of |P| in the bracket.
    let candidate = bsgs_find_multiple::<P, L>(point, curve, lo, w);

    //  Reduce to the exact order by factoring and peeling.
    let factors = trial_factor::<L>(&candidate);
    order_from_group_order::<P, L>(point, curve, &candidate, &factors)
}

// ---------------------------------------------------------------------------
// BSGS                                       (Case-B helper)
// ---------------------------------------------------------------------------
//
//  Given `P` and a bracket `[lo, lo + w²]` claimed to contain `|P|`,
//  find a value  c ∈ [lo, lo + w²]  with  [c] P = O.
//
//  Standard BSGS:
//
//    baby[j] = [j] P                 for  j = 0, 1, …, w
//    giant_i = [lo + i·w] P          for  i = 0, 1, …, w
//
//  On a collision  giant_i == baby[j]  we have  [lo + i·w - j] P = O,
//  so `c = lo + i·w - j`  is a valid multiple of the order within the
//  bracket (Phase 3 reduces it to the exact order).

fn bsgs_find_multiple<P, const L: usize>(
    point: &P,
    curve: &P::Curve,
    lo: &Uint<L>,
    w: u64,
) -> Uint<L>
where
    P: PointAdd + Eq + Hash,
{
    //  ---- Baby steps:   j · P   for  j ∈ [0, w] ----
    //
    //  When the true order |P| is smaller than w, the baby table contains
    //  duplicate keys (e.g. [0]P = [|P|]P = [2|P|]P = ...).  We keep the
    //  smallest `j` on collision via `.or_insert`; a `HashMap::insert` that
    //  overwrote would also be correct, but produces slightly larger `c`
    //  values for Phase 3 to factor.
    //
    //  Including  j = w  is needed so that the range of representable
    //  candidates  `lo + i·w − j`  covers every integer in  [lo, lo + w²].

    let mut table: HashMap<P, u64> = HashMap::with_capacity(w as usize + 1);
    let mut jp = P::identity(curve);
    for j in 0..=w {
        table.entry(jp.clone()).or_insert(j);
        jp = jp.add(point, curve);
    }

    //  ---- Giant steps: start at [lo] P, stride by [w] P ----
    //
    //  On a collision we compute  c = lo + i·w − j.  There are a few
    //  cases to handle carefully:
    //
    //   (1) `lo + i·w < j`: the (signed) candidate is negative — this
    //       just means the multiple of |P| lying in the search lattice
    //       is below `lo`.  Skip and keep walking.
    //   (2) `c = 0`: the collision witnesses a trivial relation like
    //       `[lo]P = [lo]P`.  No information about the order.  Skip.
    //   (3) `c > 0`: a valid positive multiple of |P| within at most
    //       `lo + w²` of `lo`.  Return it; Phase 3 will reduce to the
    //       exact order by trial-factor + peel.
    //
    //  Case (2) is common when |P| is small enough that the baby table
    //  contains duplicates: multiple `j` map to the same point.  The
    //  `.or_insert` above keeps the smallest `j`, which then produces
    //  `c = 0` at `i = 0`; we simply keep walking until `i` is large
    //  enough that `lo + i·w` is a genuine multiple of |P|.

    let w_uint = Uint::<L>::from(w);
    let stride = point.scalar_mul(&[w], curve);
    let mut cur = point.scalar_mul(lo.as_words(), curve);

    for i in 0..=w {
        if let Some(&j) = table.get(&cur) {
            //  Candidate c = lo + i·w − j .
            let i_w = Uint::<L>::from(i).wrapping_mul(&w_uint);
            let plus = lo.wrapping_add(&i_w);
            let j_uint = Uint::<L>::from(j);

            //  Cases (1) and (2): c ≤ 0 .  Skip.
            if j_uint < plus {
                let c = plus.wrapping_sub(&j_uint);
                if c != Uint::<L>::ZERO {
                    //  Case (3): valid positive multiple of |P|.
                    return c;
                }
            }
            //  Fall through — collision did not yield a usable c.
        }
        cur = cur.add(&stride, curve);
    }

    panic!("BSGS exhausted with no collision — bracket [lo, hi] is incorrect");
}

// ---------------------------------------------------------------------------
// Follow-up: Sutherland (2007, §5) smoothness optimisation
// ---------------------------------------------------------------------------
//
//  The present `order_in_interval` does plain BSGS over the whole
//  Hasse width and achieves  O(q^{1/4})  in the worst case.
//
//  Sutherland shows that by handling small-prime factors of  |P|
//  separately one can reduce the BSGS width to  O(\sqrt\ell_max(|P|)) ,
//  which is much less than O(q^{1/4}) whenever |P| has any small
//  prime factor. The algorithm splits |P| = m t with \ell_max(m) \leq B
//  and searches for m and t independently.
//
//  A correct implementation requires: (a) a divisibility oracle that
//  does not assume |P| divides any fixed round number, (b) careful
//  lifting between the \ell-torsion and the full group, and (c) the
//  orchestration described in Sutherland's Algorithm 5.1.
//
//  Tracked here so the TODO is discoverable from the module root.
//
//  TODO(sutherland-phase-1):
//    - Implement  order_in_interval_smooth<P, L>(...)
//    - Match the I(sqrt(\elll)_max(|P|)) bound
//    - Add reference tests against known-order curves (e.g. P-256, where
//      |E(F_p)| is prime so the new algorithm degenerates to the current
//      one and must agree bit-for-bit).
