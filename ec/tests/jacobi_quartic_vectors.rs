//! Reusable test suite for the affine Jacobi quartic implementation.
//!
//! This file is meant to live under `tests/` and be included from one small
//! field-specific wrapper. The wrapper only has to provide:
//!
//! - the ec crate path,
//! - a concrete field type implementing `FieldOps`,
//! - a lifting function `fn(i64) -> F`.
//!
//! Example:
//!
//! ```ignore
//! #[path = "jacobi_quartic_vectors.rs"]
//! mod jacobi_quartic_vectors;
//!
//! use my_field_adapter::{f, F};
//!
//! jacobi_quartic_vectors::jacobi_quartic_test_suite!(ec, F, f);
//! ```

#[macro_export]
macro_rules! jacobi_quartic_test_suite {
    ($ec:ident, $field_ty:ty, $lift:path)  => {
        use $ec::curve_jacobi_quartic::JacobiQuarticCurve;
        use $ec::point_jacobi_quartic::JacobiQuarticPoint;

        fn jqf(n: i64) -> $field_ty {
            $lift(n)
        }

        fn jq_curve() -> JacobiQuarticCurve<$field_ty> {
            // Sage fixture: GF(101), a = -1/2 = 50, d = 2 (nonsquare).
            JacobiQuarticCurve::new(jqf(50), jqf(2))
        }

        fn jq_id() -> JacobiQuarticPoint<$field_ty> {
            JacobiQuarticPoint::new(jqf(0), jqf(1))
        }

        fn jq_p() -> JacobiQuarticPoint<$field_ty> {
            JacobiQuarticPoint::new(jqf(8), jqf(94))
        }

        fn jq_q() -> JacobiQuarticPoint<$field_ty> {
            JacobiQuarticPoint::new(jqf(9), jqf(35))
        }

        #[test]
        fn jacobi_quartic_contains_known_points() {
            let curve = jq_curve();
            let id = jq_id();
            let p = jq_p();
            let q = jq_q();
            let minus_p = p.negate(&curve);

            assert!(curve.contains(&id.x, &id.y));
            assert!(curve.contains(&p.x, &p.y));
            assert!(curve.contains(&q.x, &q.y));
            assert!(curve.contains(&minus_p.x, &minus_p.y));
        }

        #[test]
        fn jacobi_quartic_identity_and_negation() {
            let curve = jq_curve();
            let id = jq_id();
            let p = jq_p();
            let minus_p = p.negate(&curve);

            assert_eq!(id, JacobiQuarticPoint::identity());
            assert_eq!(p.add(&id, &curve), p);
            assert_eq!(id.add(&p, &curve), p);
            assert_eq!(p.add(&minus_p, &curve), id);
            assert_eq!(minus_p.add(&p, &curve), id);
        }

        #[test]
        fn jacobi_quartic_known_doubling_vectors() {
            let curve = jq_curve();
            let p = jq_p();
            let q = jq_q();

            assert_eq!(p.double(&curve), JacobiQuarticPoint::new(jqf(92), jqf(35)));
            assert_eq!(q.double(&curve), JacobiQuarticPoint::new(jqf(70), jqf(99)));
        }

        #[test]
        fn jacobi_quartic_known_addition_vectors() {
            let curve = jq_curve();
            let p = jq_p();
            let q = jq_q();

            let expected = JacobiQuarticPoint::new(jqf(93), jqf(94));
            assert_eq!(p.add(&q, &curve), expected);
            assert_eq!(q.add(&p, &curve), expected);
        }

        #[test]
        fn jacobi_quartic_scalar_mul_vectors() {
            let curve = jq_curve();
            let p = jq_p();

            assert_eq!(p.scalar_mul(&[0], &curve), JacobiQuarticPoint::new(jqf(0), jqf(1)));
            assert_eq!(p.scalar_mul(&[1], &curve), JacobiQuarticPoint::new(jqf(8), jqf(94)));
            assert_eq!(p.scalar_mul(&[2], &curve), JacobiQuarticPoint::new(jqf(92), jqf(35)));
            assert_eq!(p.scalar_mul(&[3], &curve), JacobiQuarticPoint::new(jqf(46), jqf(37)));
            assert_eq!(p.scalar_mul(&[5], &curve), JacobiQuarticPoint::new(jqf(89), jqf(11)));
            assert_eq!(p.scalar_mul(&[8], &curve), JacobiQuarticPoint::new(jqf(21), jqf(50)));
            assert_eq!(p.scalar_mul(&[16], &curve), JacobiQuarticPoint::new(jqf(93), jqf(94)));
            assert_eq!(p.scalar_mul(&[17], &curve), jq_id());
        }

        #[test]
        fn jacobi_quartic_double_matches_add_self() {
            let curve = jq_curve();
            let p = jq_p();
            assert_eq!(p.double(&curve), p.add(&p, &curve));
        }

        #[test]
        fn jacobi_quartic_small_associativity_sample() {
            let curve = jq_curve();
            let p = jq_p();
            let q = jq_q();
            let r = p.double(&curve);

            assert_eq!(p.add(&q, &curve).add(&r, &curve), p.add(&q.add(&r, &curve), &curve));
        }

        #[test]
        fn jacobi_quartic_j_invariant_matches_fixture() {
            let curve = jq_curve();
            // Over GF(101): j = 59 for a = 50, d = 2.
            assert_eq!(curve.j_invariant(), jqf(59));
        }
    };
}
