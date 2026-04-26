//! Reusable test suite for the affine Jacobi intersection implementation.
//!
//! The wrapper pattern is the same as for `jacobi_quartic_vectors.rs`.

#[macro_export]
macro_rules! jacobi_intersection_test_suite {
    ($ec:ident, $field_ty:ty, $lift:path) => {
        use $ec::curve_jacobi_intersection::JacobiIntersectionCurve;
        use $ec::point_jacobi_intersection::JacobiIntersectionPoint;

        fn jif(n: i64) -> $field_ty {
            $lift(n)
        }

        fn ji_curve() -> JacobiIntersectionCurve<$field_ty> {
            // Sage fixture: GF(101), a = 2.
            JacobiIntersectionCurve::new(jif(2))
        }

        fn ji_id() -> JacobiIntersectionPoint<$field_ty> {
            JacobiIntersectionPoint::new(jif(0), jif(1), jif(1))
        }

        fn ji_p() -> JacobiIntersectionPoint<$field_ty> {
            JacobiIntersectionPoint::new(jif(5), jif(28), jif(31))
        }

        fn ji_q() -> JacobiIntersectionPoint<$field_ty> {
            JacobiIntersectionPoint::new(jif(20), jif(56), jif(98))
        }

        #[test]
        fn jacobi_intersection_contains_known_points() {
            let curve = ji_curve();
            let id = ji_id();
            let p = ji_p();
            let q = ji_q();
            let minus_p = p.negate(&curve);

            assert!(curve.contains(&id.s, &id.c, &id.d));
            assert!(curve.contains(&p.s, &p.c, &p.d));
            assert!(curve.contains(&q.s, &q.c, &q.d));
            assert!(curve.contains(&minus_p.s, &minus_p.c, &minus_p.d));
        }

        #[test]
        fn jacobi_intersection_identity_and_negation() {
            let curve = ji_curve();
            let id = ji_id();
            let p = ji_p();
            let minus_p = p.negate(&curve);

            assert_eq!(id, JacobiIntersectionPoint::identity());
            assert_eq!(p.add(&id, &curve), p);
            assert_eq!(id.add(&p, &curve), p);
            assert_eq!(p.add(&minus_p, &curve), id);
            assert_eq!(minus_p.add(&p, &curve), id);
        }

        #[test]
        fn jacobi_intersection_known_doubling_vectors() {
            let curve = ji_curve();
            let p = ji_p();
            let q = ji_q();

            assert_eq!(
                p.double(&curve),
                JacobiIntersectionPoint::new(jif(22), jif(74), jif(89))
            );
            assert_eq!(q.double(&curve), ji_p());
        }

        #[test]
        fn jacobi_intersection_known_addition_vectors() {
            let curve = ji_curve();
            let p = ji_p();
            let q = ji_q();

            let expected = JacobiIntersectionPoint::new(jif(44), jif(40), jif(88));
            assert_eq!(p.add(&q, &curve), expected);
            assert_eq!(q.add(&p, &curve), expected);
        }

        #[test]
        fn jacobi_intersection_scalar_mul_vectors() {
            let curve = ji_curve();
            let p = ji_p();

            assert_eq!(p.scalar_mul(&[0], &curve), ji_id());
            assert_eq!(p.scalar_mul(&[1], &curve), ji_p());
            assert_eq!(
                p.scalar_mul(&[2], &curve),
                JacobiIntersectionPoint::new(jif(22), jif(74), jif(89))
            );
            assert_eq!(
                p.scalar_mul(&[3], &curve),
                JacobiIntersectionPoint::new(jif(52), jif(78), jif(42))
            );
            assert_eq!(
                p.scalar_mul(&[5], &curve),
                JacobiIntersectionPoint::new(jif(57), jif(40), jif(88))
            );
            assert_eq!(
                p.scalar_mul(&[8], &curve),
                JacobiIntersectionPoint::new(jif(44), jif(40), jif(88))
            );
            assert_eq!(p.scalar_mul(&[13], &curve), ji_id());
        }

        #[test]
        fn jacobi_intersection_double_matches_add_self() {
            let curve = ji_curve();
            let p = ji_p();
            assert_eq!(p.double(&curve), p.add(&p, &curve));
        }

        #[test]
        fn jacobi_intersection_small_associativity_sample() {
            let curve = ji_curve();
            let p = ji_p();
            let q = ji_q();
            let r = p.double(&curve);

            assert_eq!(
                p.add(&q, &curve).add(&r, &curve),
                p.add(&q.add(&r, &curve), &curve)
            );
        }
    };
}
