//! Integration tests for `ec::point_order_ops`.
//!
//! Covers Sutherland (2007):
//!
//!   * [`order_from_group_order`]  — fast path when `N = #E(F_q)` is known
//!   * [`order_in_interval`]       — generic BSGS on the Hasse interval
//!
//! # Fixture
//!
//!   * Two distinct primes (`2` and `3`) with a squared factor → exercises
//!     the inner peeling loop of [`order_from_group_order`].
//!   * Full order spectrum represented in the group:
//!
//!     | order | count | role                 |
//!     |-------|-------|----------------------|
//!     | 1     | 1     | identity `O`         |
//!     | 2     | 1     | 2-torsion            |
//!     | 3     | 2     | 3-torsion            |
//!     | 6     | 2     | 6-torsion            |
//!     | 9     | 6     | inside Hasse bracket |
//!     | 18    | 6     | generators           |
//!
//!   * Hasse interval `[q+1-⌈2√q⌉, q+1+⌈2√q⌉] = [10, 26]` for `q = 17`,
//!     which contains orders `18` and `9` — so BSGS (`w = ⌈√16⌉ = 4`) has
//!     real work to do.
//!
//! Every test that needs the "true" order of a point gets it via brute-force
//! repeated addition (legal here because |P| ≤ 18), so our reference is
//! independent of the algorithms under test.

use crypto_bigint::{Uint, const_prime_monty_params};

use fp::fp_element::FpElement;

use ec::curve_edwards::EdwardsCurve;
use ec::curve_weierstrass::WeierstrassCurve;
use ec::num_utils::{SmallPrimes, product_of_factors, trial_factor};
use ec::point_edwards::EdwardsPoint;
use ec::point_order_ops::{order_from_group_order, order_in_interval};
use ec::point_weierstrass::AffinePoint;

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------

const_prime_monty_params!(Fp17Mod, Uint<1>, "0000000000000011", 3);
type F17 = FpElement<Fp17Mod, 1>;

fn fp(n: u64) -> F17 {
    F17::from_u64(n)
}

/// `y² = x³ + 1` over `F₁₇`
fn curve() -> WeierstrassCurve<F17> {
    WeierstrassCurve::new_short(fp(0), fp(1))
}

/// Expected group order `#E(F₁₇) = 18`
fn group_order() -> Uint<1> {
    Uint::<1>::from(18u64)
}

/// Expected factorization `18 = 2 · 3²`
fn group_order_factors() -> Vec<(Uint<1>, u32)> {
    vec![(Uint::<1>::from(2u64), 1), (Uint::<1>::from(3u64), 2)]
}

/// Hasse bracket for `q = 17`:  `[10, 26]` (width 16, `w = 4`).
fn hasse_interval() -> (Uint<1>, Uint<1>) {
    (Uint::<1>::from(10u64), Uint::<1>::from(26u64))
}

/// All affine points on the fixture curve (brute force over F₁₇).
///
/// Does **not** include the point at infinity — callers append it if
/// they want to test the identity too.
fn all_affine_points() -> Vec<(u64, u64)> {
    let p = 17u64;
    let mut pts = Vec::new();
    for x in 0..p {
        //  rhs = x³ + 1  (mod p)
        let rhs = (x.wrapping_mul(x) % p * x % p + 1) % p;
        for y in 0..p {
            if (y.wrapping_mul(y)) % p == rhs {
                pts.push((x, y));
            }
        }
    }
    pts
}

/// All points including `O`, as `AffinePoint<F₁₇>` values.
fn all_points() -> Vec<AffinePoint<F17>> {
    let mut v: Vec<AffinePoint<F17>> = all_affine_points()
        .into_iter()
        .map(|(x, y)| AffinePoint::new(fp(x), fp(y)))
        .collect();
    v.push(AffinePoint::<F17>::identity());
    v
}

/// Brute-force order: smallest `k ≥ 1` with `[k] P = O`, via repeated
/// addition.  Only acceptable because |P| ≤ 18 on this fixture.
fn brute_force_order(p: &AffinePoint<F17>, c: &WeierstrassCurve<F17>) -> u64 {
    if p.is_identity() {
        return 1;
    }
    let mut acc = *p;
    let mut k: u64 = 1;
    while !acc.is_identity() {
        acc = acc.add(p, c);
        k += 1;
        assert!(k < 100, "runaway in brute-force order");
    }
    k
}

// ---------------------------------------------------------------------------
// Fixture self-tests — fail early if the curve equation ever changes
// ---------------------------------------------------------------------------

#[test]
fn fixture_has_expected_point_count() {
    //  17 affine points on y² = x³ + 1 over F₁₇, plus the point at
    //  infinity, for a total of #E(F₁₇) = 18.
    let n = all_affine_points().len() as u64 + 1;
    assert_eq!(
        n, 18,
        "fixture curve should have exactly 18 points, got {}",
        n
    );
}

#[test]
fn fixture_has_expected_order_distribution() {
    use std::collections::BTreeMap;
    let c = curve();
    let mut dist: BTreeMap<u64, u32> = BTreeMap::new();
    for p in all_points() {
        *dist.entry(brute_force_order(&p, &c)).or_insert(0) += 1;
    }
    //  Expected: {1:1, 2:1, 3:2, 6:2, 9:6, 18:6}
    let expected: BTreeMap<u64, u32> = [(1, 1), (2, 1), (3, 2), (6, 2), (9, 6), (18, 6)]
        .into_iter()
        .collect();
    assert_eq!(dist, expected, "order distribution regression");
}

// ---------------------------------------------------------------------------
// num_utils at the integration level
// ---------------------------------------------------------------------------

#[test]
fn num_utils_sieve_agrees_with_hand_list() {
    let primes_up_to_30: Vec<u64> = SmallPrimes::up_to(30).collect();
    assert_eq!(primes_up_to_30, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
}

#[test]
fn num_utils_factors_group_order() {
    let n = Uint::<1>::from(18u64);
    let f = trial_factor(&n);
    assert_eq!(f, group_order_factors());
    assert_eq!(product_of_factors(&f), n);
}

#[test]
fn num_utils_round_trips_every_small_int() {
    //  For every n in [1, 50], trial_factor then product_of_factors must
    //  return n.  Catches off-by-ones in either direction.
    for n_u64 in 1u64..=50 {
        let n = Uint::<1>::from(n_u64);
        let f = trial_factor(&n);
        assert_eq!(
            product_of_factors(&f),
            n,
            "round-trip failed on n = {}",
            n_u64
        );
    }
}

// ---------------------------------------------------------------------------
// Case A:  order_from_group_order
// ---------------------------------------------------------------------------

#[test]
fn case_a_identity_has_order_one() {
    let c = curve();
    let n = group_order();
    let factors = group_order_factors();

    let id = AffinePoint::<F17>::identity();
    let ord = order_from_group_order(&id, &c, &n, &factors);
    assert_eq!(ord, Uint::<1>::from(1u64), "order of O should be 1");
}

#[test]
fn case_a_every_point_matches_brute_force() {
    let c = curve();
    let n = group_order();
    let factors = group_order_factors();

    for p in all_points() {
        let expected = brute_force_order(&p, &c);
        let got = order_from_group_order(&p, &c, &n, &factors);
        assert_eq!(
            got,
            Uint::<1>::from(expected),
            "point {}: brute-force = {}, order_from_group_order = {:?}",
            p,
            expected,
            got,
        );
        //  Sanity: every order must divide 18.
        assert_eq!(18 % expected, 0, "order {} does not divide 18", expected);
    }
}

#[test]
fn case_a_finds_a_generator() {
    //  At least one point must be reported as order 18 (i.e. a generator).
    let c = curve();
    let n = group_order();
    let factors = group_order_factors();

    let mut found = false;
    for p in all_points() {
        let ord = order_from_group_order(&p, &c, &n, &factors);
        if ord == Uint::<1>::from(18u64) {
            found = true;
            break;
        }
    }
    assert!(found, "expected at least one generator of order 18");
}

const_prime_monty_params!(
    Fp25519Mod,
    Uint<4>,
    "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed",
    2
);
type F25519 = FpElement<Fp25519Mod, 4>;

fn fp25519(words: [u64; 4]) -> F25519 {
    F25519::from_words(words)
}

/// Untwisted Edwards model birationally equivalent to Curve25519:
///
///     x² + y² = 1 + d x² y²,   d = (A - 2)/(A + 2),   A = 486662.
///
/// This is the odd-characteristic Edwards model implemented by
/// `ec::curve_edwards::EdwardsCurve`.
fn curve25519_edwards() -> EdwardsCurve<F25519> {
    // d = 121665 / 121666 mod p
    EdwardsCurve::new(fp25519([
        9949773421441615690,
        18415207549394364244,
        8302596497594521447,
        3313685130627776908,
    ]))
}

/// Prime-order subgroup size of the Curve25519 base point.
fn curve25519_subgroup_order() -> Uint<4> {
    Uint::<4>::from_words([
        6346243789798364141,
        1503914060200516822,
        0,
        1152921504606846976,
    ])
}

/// Full Curve25519 group order = 8 · ℓ.
fn curve25519_group_order() -> Uint<4> {
    Uint::<4>::from_words([
        13876462170967809896,
        12031312481604134578,
        0,
        9223372036854775808,
    ])
}

/// Factorization of `#E(F_p) = 2³ · ℓ`, with `ℓ` prime.
fn curve25519_group_order_factors() -> Vec<(Uint<4>, u32)> {
    vec![
        (Uint::<4>::from(2u64), 3),
        (curve25519_subgroup_order(), 1),
    ]
}

/// Curve25519 Montgomery base point `u = 9`, mapped into the untwisted
/// Edwards model above.
///
/// We hard-code one of the two sign choices for `x`; both have the same
/// order.  The `y` coordinate is `(u - 1)/(u + 1) = 4/5`.
fn curve25519_basepoint_edwards() -> EdwardsPoint<F25519> {
    EdwardsPoint::new(
        fp25519([
            14735460775765454235,
            6707206447363776329,
            440012832734768475,
            1556671287782030463,
        ]),
        fp25519([
            7378697629483820632,
            7378697629483820646,
            7378697629483820646,
            7378697629483820646,
        ]),
    )
}

#[test]
fn curve25519_fixture_is_consistent() {
    let c = curve25519_edwards();
    let p = curve25519_basepoint_edwards();
    let ell = curve25519_subgroup_order();

    assert!(
        c.contains(&p.x, &p.y),
        "mapped Curve25519 base point must lie on the Edwards model",
    );

    let killed = p.scalar_mul(ell.as_words(), &c);
    assert!(
        killed.is_identity(),
        "[ℓ]P must be O for the Curve25519 base point",
    );

    let small = p.scalar_mul(&[8u64], &c);
    assert!(
        !small.is_identity(),
        "base point should not lie in the small cofactor subgroup",
    );
}

#[test]
fn curve25519_case_a_recovers_prime_subgroup_order() {
    let c = curve25519_edwards();
    let p = curve25519_basepoint_edwards();
    let n = curve25519_group_order();
    let factors = curve25519_group_order_factors();
    let ell = curve25519_subgroup_order();

    let ord = order_from_group_order(&p, &c, &n, &factors);
    assert_eq!(
        ord, ell,
        "Curve25519 base point should have prime order ℓ inside #E(F_p)=8·ℓ",
    );
}

#[test]
fn curve25519_case_a_order_is_minimal() {
    let c = curve25519_edwards();
    let p = curve25519_basepoint_edwards();
    let n = curve25519_group_order();
    let factors = curve25519_group_order_factors();
    let ell = curve25519_subgroup_order();

    let ord = order_from_group_order(&p, &c, &n, &factors);
    assert_eq!(ord, ell);

    let half = ord.wrapping_shr_vartime(1);
    let maybe_killed = p.scalar_mul(half.as_words(), &c);
    assert!(
        !maybe_killed.is_identity(),
        "[ℓ/2]P must not be O when |P| = ℓ is an odd prime",
    );
}

// ---------------------------------------------------------------------------
// Case B:  order_in_interval
// ---------------------------------------------------------------------------

#[test]
fn case_b_matches_brute_force_for_points_in_hasse() {
    //  For every point whose order lies in the Hasse interval, BSGS must
    //  recover it exactly.  Orders on this fixture inside [10, 26] are
    //  {18, 9} — that's 12 out of 18 points.
    let c = curve();
    let (lo, hi) = hasse_interval();

    let mut tested = 0;
    for p in all_points() {
        let true_order = brute_force_order(&p, &c);
        if !(10u64..=26).contains(&true_order) {
            continue;
        }
        tested += 1;

        let got = order_in_interval(&p, &c, &lo, &hi);
        assert_eq!(
            got,
            Uint::<1>::from(true_order),
            "order_in_interval: point {}: expected {}, got {:?}",
            p,
            true_order,
            got,
        );
    }
    assert!(
        tested >= 6,
        "expected ≥6 points in Hasse, tested {}",
        tested
    );
}

#[test]
fn case_b_agrees_with_case_a_on_wide_bracket() {
    //  With a bracket wide enough to cover every possible order
    //  ([1, 18]), Case A and Case B must agree on every point.
    let c = curve();
    let n = group_order();
    let factors = group_order_factors();
    let lo = Uint::<1>::from(1u64);
    let hi = Uint::<1>::from(18u64);

    for p in all_points() {
        let ord_a = order_from_group_order(&p, &c, &n, &factors);
        let ord_b = order_in_interval(&p, &c, &lo, &hi);
        assert_eq!(
            ord_a, ord_b,
            "Case A ≠ Case B for point {}: A={:?}, B={:?}",
            p, ord_a, ord_b,
        );
    }
}

#[test]
fn case_b_tight_brackets_each_order() {
    //  For each achievable order `k` on this curve, check that Case B
    //  returns `k` when called with the tight bracket [k, k].
    //  A bracket of width 0 gives w = 0, so this also exercises the
    //  degenerate (but legal) BSGS boundary.
    let c = curve();
    for k in [1u64, 2, 3, 6, 9, 18] {
        let lo = Uint::<1>::from(k);
        let hi = Uint::<1>::from(k);
        for p in all_points() {
            if brute_force_order(&p, &c) != k {
                continue;
            }
            let got = order_in_interval(&p, &c, &lo, &hi);
            assert_eq!(
                got,
                Uint::<1>::from(k),
                "tight bracket [{k},{k}] on order-{k} point {p}: got {got:?}",
            );
            //  One sample per order is enough.
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// Post-conditions (property-style)
// ---------------------------------------------------------------------------

#[test]
fn returned_order_annihilates_point() {
    //  Whatever order we return, [order] P must be O.
    //  Catches off-by-one bugs in the peeling loop.
    let c = curve();
    let n = group_order();
    let factors = group_order_factors();

    for p in all_points() {
        let ord = order_from_group_order(&p, &c, &n, &factors);
        let killed = p.scalar_mul(ord.as_words(), &c);
        assert!(
            killed.is_identity(),
            "[order] P must be O, but wasn't for {p}, order={ord:?}",
        );
    }
}

#[test]
fn returned_order_is_minimal() {
    //  For every prime `ℓ | order(P)`, [order/ℓ] P must NOT be O.
    //  This distinguishes "exact order" from "a multiple of the order".
    let c = curve();
    let n = group_order();
    let factors = group_order_factors();

    for p in all_points() {
        let ord = order_from_group_order(&p, &c, &n, &factors);
        let ord_u64 = ord.as_words()[0];
        if ord_u64 <= 1 {
            continue;
        }
        for (prime, _) in &factors {
            let prime_u64 = prime.as_words()[0];
            if ord_u64 % prime_u64 != 0 {
                continue;
            }
            let divisor = ord_u64 / prime_u64;
            let maybe_killed = p.scalar_mul(&[divisor], &c);
            assert!(
                !maybe_killed.is_identity(),
                "[{divisor}] P killed point {p}, so claimed order {ord_u64} is not minimal",
            );
        }
    }
}

#[test]
fn case_b_returned_order_annihilates_point() {
    //  Same post-condition, but for Case B on the Hasse bracket.
    let c = curve();
    let (lo, hi) = hasse_interval();

    for p in all_points() {
        if !(10u64..=26).contains(&brute_force_order(&p, &c)) {
            continue;
        }
        let ord = order_in_interval(&p, &c, &lo, &hi);
        let killed = p.scalar_mul(ord.as_words(), &c);
        assert!(
            killed.is_identity(),
            "Case B [order] P must be O, but wasn't for {p}, order={ord:?}",
        );
    }
}
