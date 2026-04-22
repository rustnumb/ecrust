//! Isogenies for montgomery curves
//!
//! AT THE MOMENT THIS WORKS ONLY IN ODD CHARACTERISTIC!!!


use subtle::{ConstantTimeEq, Choice};

use fp::field_ops::FieldOps;
use ec::curve_montgomery::MontgomeryCurve;
use ec::point_montgomery::KummerPoint;
use fp::{ref_field_fns, ref_field_impl, ref_field_trait_impl};

use primal::is_prime;

#[derive(Debug, Clone)]
pub struct PointAndDegree<F: FieldOps + Copy> {
    pub pt: KummerPoint<F>,
    pub degree: usize,
    pub degree_known: Choice
}

#[derive(Debug, Clone)]
pub struct MontgomeryIsogeny<F: FieldOps + Copy> {
    pub domain: MontgomeryCurve<F>,
    pub codomain: MontgomeryCurve<F>,
    // x-coordinate of the generator of the kernel K = <G> + degree of the isogeny
    pub pt_deg: PointAndDegree<F>
}


// ---------------------------------------------------------------------------
// Manual trait impls
// ---------------------------------------------------------------------------

ref_field_trait_impl! {
    impl<F> PartialEq for PointAndDegree<F> {
        fn eq(&self, other: &Self) -> bool {
            self.pt == other.pt
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + ConstantTimeEq> Eq for PointAndDegree<F> { }
}

ref_field_trait_impl! {
    impl<F> PartialEq for MontgomeryIsogeny<F> {
        fn eq(&self, other: &Self) -> bool {
            self.domain == other.domain && self.codomain == other.codomain && self.pt_deg == other.pt_deg
        }
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + ConstantTimeEq> Eq for MontgomeryIsogeny<F> { }
}


// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

ref_field_impl! {
    impl<F> PointAndDegree<F> {
        pub fn new_known_degree(pt: KummerPoint<F>, degree: usize) -> Self {
            Self{pt, degree, degree_known: Choice::from(1)}
        }

        pub fn new_unknown_degree(pt: KummerPoint<F>) -> Self {
            Self{pt, degree: 0, degree_known: Choice::from(0)}
        }
    }
}


ref_field_impl! {
    impl<F> MontgomeryIsogeny<F> {
        pub fn new(domain: MontgomeryCurve<F>, pt_deg: PointAndDegree<F>) -> Self {
            let x_pts = kernel_rge_prime_degree(&domain, &pt_deg);
            let (a_codomain, b_codomain) = new_codomain(&domain, &x_pts);
            Self{
                domain,
                codomain: MontgomeryCurve::new(a_codomain, b_codomain),
                pt_deg
            }
        }

        pub fn degree(domain: &MontgomeryCurve<F>, pt: KummerPoint<F>) -> usize {
            todo!()
        }

    }
}


// ---------------------------------------------------------------------------
// Utilities and helper functions
// ---------------------------------------------------------------------------

ref_field_fns! {
    fn criss_cross<F>(a: &F, b: &F, c: &F, d: &F) -> (F, F) {
        let t1 = a * d;
        let t2 = b * c;
        (&t1 + &t2, &t1 - &t2)
    }

        // Costello-Hisil algorithm
    fn kernel_rge_prime_degree<F>(domain: &MontgomeryCurve<F>, pt_deg: &PointAndDegree<F>) -> Vec<KummerPoint<F>> {
        assert!(bool::from(pt_deg.degree_known), "Degree needs to be known!");
        let l = pt_deg.degree as u64;
        assert!(is_prime(l) & (l != 2));

        let mut res = Vec::new();
        let pt = pt_deg.pt;

        res.push(pt.clone());
        res.push(pt.xdouble(domain));

        for i in 3..=(l-1)/2 {
            let new_pt = res[(i-1) as usize].xadd(&pt, &res[(i-2) as usize]);
            res.push(new_pt);
        }

        res
    }

    fn new_codomain<F>(domain: &MontgomeryCurve<F>, x_pts: &Vec<KummerPoint<F>>) -> (F, F) {
        let d = x_pts.len();

        let mut sigma = F::zero();
        let mut sigma_tilde = F::zero();
        let mut pi = F::one();

        for i in 0..d {
            let xi = x_pts[i].to_x().unwrap();
            sigma = &sigma + &xi;
            sigma_tilde = &sigma_tilde + &<F as FieldOps>::invert(&xi).unwrap();
            pi = &pi * &xi;
        }

        let diff_sigmas = &sigma - &sigma_tilde;
        let six_diff_sigmas = &F::from_u64(6) * &diff_sigmas;
        let tmp = &six_diff_sigmas + &domain.a;
        let pi_sq = <F as FieldOps>::square(&pi);

        let new_a = &tmp * &pi_sq;
        let new_b = &domain.b * &pi_sq;

        (new_a, new_b)
    }

    fn evaluate<F>(x_pts: &Vec<KummerPoint<F>>, eval_pt: KummerPoint<F>) -> KummerPoint<F> {
        let mut new_x = eval_pt.x.clone();
        let mut new_z = eval_pt.z.clone();

        for pt in x_pts {
            let tmp_x = &(&eval_pt.x * &pt.x) - &(&eval_pt.z * &pt.z);
            let tmp_z = &(&eval_pt.x * &pt.z) - &(&eval_pt.z * &pt.x);

            new_x = &new_x * &<F as FieldOps>::square(&tmp_x);
            new_z = &new_z * &<F as FieldOps>::square(&tmp_z);
        }

        KummerPoint{x: new_x, z: new_z}
    }
}


// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

