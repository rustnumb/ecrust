//! Small-integer utilities shared by the `ec` crate.
//!
//! These are deliberately simple, `alloc`-using helpers intended for
//! auxiliary computations like order determination and point-count
//! post-processing — **not** for the hot inner loops of scalar
//! multiplication. They operate on `u64` for compactness and on
//! [`Uint<L>`](crypto_bigint::Uint) for interop with the rest of the
//! stack.
//!
//! # Contents
//!
//! | Function                | What it does                                       |
//! |-------------------------|----------------------------------------------------|
//! | [`SmallPrimes::up_to`]  | Sieve of Eratosthenes, lazy iterator over primes   |
//! | [`trial_factor`]        | Prime factorization by trial division              |
//! | [`product_of_factors`]  | Multiplies out `[(pᵢ, eᵢ)]` (debug-assert helper)  |

use crypto_bigint::{NonZero, Uint};

// ---------------------------------------------------------------------------
// SmallPrimes — sieve of Eratosthenes over [2, bound]
// ---------------------------------------------------------------------------

/// Iterator yielding the primes `ℓ ≤ bound`, in increasing order.
///
/// Construction is $O(B \log \log B)$ time and $O(B / 8)$ memory (one
/// bit per odd integer).  The iterator itself then yields each prime
/// in amortised $O(1)$ by advancing a cursor over the sieve.
///
/// Intended for `bound` up to roughly $2^{20}$; beyond that you should
/// prefer a segmented sieve.
///
/// # Example
///
/// ```text
/// let small: Vec<u64> = SmallPrimes::up_to(20).collect();
/// // small == [2, 3, 5, 7, 11, 13, 17, 19]
/// ```
pub struct SmallPrimes {
    //  `sieve[i] = true`  ⇔  (2i + 1)  is composite,  except for the
    //  explicit special case of 2 handled by `yielded_two`.
    //  We sieve only odd numbers to halve memory.
    sieve: Vec<bool>,
    //  Cursor in odd-number space: next candidate is  2·cursor + 1.
    cursor: usize,
    //  Inclusive upper bound in odd-number space.
    odd_limit: usize,
    //  Whether we have already emitted the prime 2.
    yielded_two: bool,
}

impl SmallPrimes {
    /// Build the sieve for primes `≤ bound`.
    pub fn up_to(bound: u64) -> Self {
        let bound_us = bound as usize;

        //  Odd-number space: index i ↔ value 2i + 1.
        //  Max odd value we need is `bound` (if bound is odd) or
        //  `bound - 1`.  The index for (2i+1 = bound) is i = (bound-1)/2.
        let odd_limit = if bound < 3 { 0 } else { (bound_us - 1) / 2 };

        let mut sieve = vec![false; odd_limit + 1];
        //  Index 0 ↔ value 1, which is not prime.
        sieve[0] = true;

        //  Cross out odd composites.  For each prime p = 2i+1 with p·p ≤ bound,
        //  mark  p·p, p·(p+2), p·(p+4), … .  In odd-index space the step is p.
        let mut i = 1usize;
        while {
            let p = 2 * i + 1;
            p.saturating_mul(p) <= bound_us
        } {
            if !sieve[i] {
                let p = 2 * i + 1;
                //  First composite to cross out is p*p; its index is (p*p - 1)/2.
                let mut j = (p * p - 1) / 2;
                while j <= odd_limit {
                    sieve[j] = true;
                    j += p; // step  2p  in value-space  =  p  in odd-index space
                }
            }
            i += 1;
        }

        Self {
            sieve,
            cursor: 0,
            odd_limit,
            yielded_two: bound < 2,
        }
    }
}

impl Iterator for SmallPrimes {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        //  Emit 2 first.
        if !self.yielded_two {
            self.yielded_two = true;
            return Some(2);
        }

        //  Then scan odd indices ≥ 1  (skipping index 0 which is the value 1).
        if self.cursor == 0 {
            self.cursor = 1;
        }
        while self.cursor <= self.odd_limit {
            if !self.sieve[self.cursor] {
                let value = (2 * self.cursor + 1) as u64;
                self.cursor += 1;
                return Some(value);
            }
            self.cursor += 1;
        }
        None
    }
}

// ---------------------------------------------------------------------------
// trial_factor — `n = ∏ pᵢ^eᵢ`  via trial division
// ---------------------------------------------------------------------------

/// Trial-factor `n` into `[(pᵢ, eᵢ)]` with `pᵢ` prime and `∏ pᵢ^eᵢ = n`.
///
/// # Algorithm
///
/// Sieve all primes up to $\sqrt n$ and divide them out one at a time.
/// Any residue greater than 1 after the loop is itself a prime factor
/// (because it has no prime factor $\le \sqrt n$).
///
/// # Complexity
///
/// $O(\pi(\sqrt n) \cdot M)$ where $M$ is the cost of a big-int division.
///
/// # When this is acceptable
///
/// Trial division is only sensible when $n$ is known to be small or
/// smooth — which is exactly the situation in order computations
/// (Sutherland's Phase 3 residue is bounded by $\sqrt{\text{Hasse width}}$).
///
/// # Returns
///
/// Empty vector for `n = 1`.  Otherwise a list of `(prime, exponent)`
/// pairs with prime factors in ascending order.
pub fn trial_factor<const L: usize>(n: &Uint<L>) -> Vec<(Uint<L>, u32)> {
    if n == &Uint::<L>::ONE {
        return Vec::new();
    }

    let mut residue = *n;
    let mut out: Vec<(Uint<L>, u32)> = Vec::new();

    //  Bound the sieve by  ⌈√n⌉ + 1, capped at u64::MAX so we can sieve
    //  with u64 primes.  If n exceeds 2^128 this becomes impractical —
    //  but for the intended uses (factor residues of size ≈ √(Hasse width))
    //  this is always fine.
    let sqrt_n = n.floor_sqrt_vartime();
    let sqrt_cap: u64 = if sqrt_n.as_words()[1..].iter().all(|&w| w == 0) {
        sqrt_n.as_words()[0].saturating_add(1)
    } else {
        u64::MAX
    };

    for p in SmallPrimes::up_to(sqrt_cap) {
        let p_uint = Uint::<L>::from(p);
        let p_nz = NonZero::new(p_uint).expect("p ≥ 2 is nonzero");

        let mut exp: u32 = 0;
        loop {
            let (q, r) = residue.div_rem(&p_nz);
            if bool::from(r.is_zero()) {
                residue = q;
                exp += 1;
            } else {
                break;
            }
        }
        if exp > 0 {
            out.push((p_uint, exp));
            if residue == Uint::<L>::ONE {
                return out;
            }
        }
    }

    //  Any residue > 1 is itself prime.
    if residue != Uint::<L>::ONE {
        out.push((residue, 1));
    }
    out
}

// ---------------------------------------------------------------------------
// product_of_factors — inverse of `trial_factor`, for consistency checks
// ---------------------------------------------------------------------------

/// Multiplies out `[(pᵢ, eᵢ)]` to produce `∏ pᵢ^eᵢ`.
///
/// Used mainly in `debug_assert!` blocks to check that a caller-supplied
/// factorization matches the claimed product.
pub fn product_of_factors<const L: usize>(factors: &[(Uint<L>, u32)]) -> Uint<L> {
    let mut acc = Uint::<L>::ONE;
    for (p, e) in factors {
        for _ in 0..*e {
            acc = acc.wrapping_mul(p);
        }
    }
    acc
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U64;

    #[test]
    fn sieve_small() {
        let got: Vec<u64> = SmallPrimes::up_to(20).collect();
        assert_eq!(got, vec![2, 3, 5, 7, 11, 13, 17, 19]);
    }

    #[test]
    fn sieve_edge_cases() {
        assert_eq!(SmallPrimes::up_to(0).collect::<Vec<_>>(), vec![]);
        assert_eq!(SmallPrimes::up_to(1).collect::<Vec<_>>(), vec![]);
        assert_eq!(SmallPrimes::up_to(2).collect::<Vec<_>>(), vec![2]);
        assert_eq!(SmallPrimes::up_to(3).collect::<Vec<_>>(), vec![2, 3]);
    }

    #[test]
    fn factor_small_composite() {
        //  360 = 2^3 · 3^2 · 5
        let n = U64::from(360u64);
        let f = trial_factor(&n);
        assert_eq!(
            f,
            vec![
                (U64::from(2u64), 3),
                (U64::from(3u64), 2),
                (U64::from(5u64), 1),
            ]
        );
        assert_eq!(product_of_factors(&f), n);
    }

    #[test]
    fn factor_prime() {
        //  101 is prime.
        let n = U64::from(101u64);
        let f = trial_factor(&n);
        assert_eq!(f, vec![(U64::from(101u64), 1)]);
    }

    #[test]
    fn factor_one() {
        let n = U64::ONE;
        let f: Vec<(U64, u32)> = trial_factor(&n);
        assert!(f.is_empty());
    }

    #[test]
    fn factor_prime_squared() {
        //  49 = 7^2 — trips the "residue > sqrt" branch if exp logic is off.
        let n = U64::from(49u64);
        let f = trial_factor(&n);
        assert_eq!(f, vec![(U64::from(7u64), 2)]);
    }
}
