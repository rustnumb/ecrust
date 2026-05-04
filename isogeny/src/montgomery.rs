//! Isogenies for montgomery curves
//!
//! AT THE MOMENT THIS WORKS ONLY IN ODD CHARACTERISTIC!!!


use subtle::{ConstantTimeEq, Choice, ConditionallySelectable};

use fp::field_ops::{FieldOps, FieldRandom};
use ec::curve_montgomery::MontgomeryCurve;
use ec::point_montgomery::KummerPoint;
use fp::{ref_field_fns, ref_field_impl, ref_field_trait_impl};

use primal::is_prime;
use crate::isogeny_ops::IsogenyOps;

#[derive(Debug, Clone, Copy)]
pub struct PointAndDegree<F: FieldOps + Copy> {
    pub pt: KummerPoint<F>,
    pub degree: u64,
    pub degree_known: Choice
}

#[derive(Debug, Clone, Copy)]
pub struct MontgomeryIsogeny<F: FieldOps + Copy> {
    pub domain: MontgomeryCurve<F>,
    pub codomain: MontgomeryCurve<F>,
    // x-coordinate of the generator of the kernel K = <G> & degree of the isogeny
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
        pub fn new_known_degree(pt: KummerPoint<F>, degree: u64) -> Self {
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
            let (pt_deg, x_pts) = if !bool::from(pt_deg.degree_known) {
                let degree = Self::degree(&domain, &pt_deg.pt);
                let new_pt_deg = PointAndDegree::new_known_degree(pt_deg.pt, degree);
                let x_pts = kernel_rge_prime_degree(&domain, &new_pt_deg);
                (new_pt_deg, x_pts)
            } else {
                let x_pts = kernel_rge_prime_degree(&domain, &pt_deg);
                (pt_deg, x_pts)
            };

            let (a_codomain, b_codomain) = new_codomain(&domain, &x_pts);

            Self{
                domain,
                codomain: MontgomeryCurve::new(a_codomain, b_codomain),
                pt_deg
            }
        }

        pub fn degree(domain: &MontgomeryCurve<F>, pt: &KummerPoint<F>) -> u64 {
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
        let d = pt_deg.degree as u64;
        assert!(is_prime(d) & (d != 2));

        let mut res = Vec::new();
        let pt = pt_deg.pt;

        res.push(pt.clone());
        res.push(pt.xdouble(domain));

        for i in 3..=(d-1)/2 {
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

        let diff_sigmas = &sigma_tilde - &sigma;
        let six_diff_sigmas = &F::from_u64(6) * &diff_sigmas;
        let tmp = &six_diff_sigmas + &domain.a;
        let pi_sq = <F as FieldOps>::square(&pi);

        let new_a = &tmp * &pi_sq;
        let new_b = &domain.b * &pi_sq;

        (new_a, new_b)
    }

    fn evaluate<F>(x_pts: &Vec<KummerPoint<F>>, eval_pt: &KummerPoint<F>) -> KummerPoint<F> {
        let mut new_x = eval_pt.x.clone();
        let mut new_z = eval_pt.z.clone();

        let u = &eval_pt.x;
        let v = &eval_pt.z;

        for pt_q in x_pts {
            let xq = &pt_q.x;
            let zq = &pt_q.z;
            let tmp_x = &(u * xq) - &(v * zq);
            let tmp_z = &(u * zq) - &(v * xq);

            new_x = &new_x * &<F as FieldOps>::square(&tmp_x);
            new_z = &new_z * &<F as FieldOps>::square(&tmp_z);
        }

        KummerPoint{x: new_x, z: new_z}
    }
}


// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------

impl<F: FieldOps + Copy> ConditionallySelectable for PointAndDegree<F> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            pt: KummerPoint::conditional_select(&a.pt, &b.pt, choice),
            degree: u64::conditional_select(&a.degree, &b.degree, choice),
            degree_known: Choice::conditional_select(&a.degree_known, &b.degree_known, choice)
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        KummerPoint::conditional_assign(&mut self.pt, &other.pt, choice);
        u64::conditional_assign(&mut self.degree, &other.degree, choice);
        Choice::conditional_assign(&mut self.degree_known, &other.degree_known, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        KummerPoint::conditional_swap(&mut a.pt, &mut b.pt, choice);
        u64::conditional_swap(&mut a.degree, &mut b.degree, choice);
        Choice::conditional_swap(&mut a.degree_known, &mut b.degree_known, choice);
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + Copy + ConstantTimeEq> ConstantTimeEq for PointAndDegree<F> {
        fn ct_eq(&self, other: &Self) -> Choice {
            self.pt.ct_eq(&other.pt)
        }

        fn ct_ne(&self, other: &Self) -> Choice {
            !self.ct_eq(other)
        }
    }
}

impl<F> ConditionallySelectable for MontgomeryIsogeny<F>
where
    F: FieldOps + Copy,
{
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            domain: MontgomeryCurve::conditional_select(&a.domain, &b.domain, choice),
            codomain: MontgomeryCurve::conditional_select(&a.codomain, &b.codomain, choice),
            pt_deg: PointAndDegree::conditional_select(&a.pt_deg, &b.pt_deg, choice),
        }
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        MontgomeryCurve::conditional_assign(&mut self.domain, &other.domain, choice);
        MontgomeryCurve::conditional_assign(&mut self.codomain, &other.codomain, choice);
        PointAndDegree::conditional_assign(&mut self.pt_deg, &other.pt_deg, choice);
    }

    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        MontgomeryCurve::conditional_swap(&mut a.domain, &mut b.domain, choice);
        MontgomeryCurve::conditional_swap(&mut a.codomain, &mut b.codomain, choice);
        PointAndDegree::conditional_swap(&mut a.pt_deg, &mut b.pt_deg, choice);
    }
}

ref_field_trait_impl! {
    impl<F: FieldOps + Copy + ConstantTimeEq> ConstantTimeEq for MontgomeryIsogeny<F> {
        fn ct_eq(&self, other: &Self) -> Choice {
            self.domain.ct_eq(&other.domain)
                & self.codomain.ct_eq(&other.codomain)
                & self.pt_deg.ct_eq(&other.pt_deg)
        }

        fn ct_ne(&self, other: &Self) -> Choice {
            !self.ct_eq(other)
        }
    }
}


// ---------------------------------------------------------------------------
// Isogeny predicates
// ---------------------------------------------------------------------------

ref_field_trait_impl! {
    impl<F: FieldOps + FieldRandom> IsogenyOps for MontgomeryIsogeny<F> {
        type BaseField = F;
        type DomainCurve = MontgomeryCurve<F>;
        type CodomainCurve = MontgomeryCurve<F>;
        type DomainPoint = KummerPoint<F>;
        type CodomainPoint = KummerPoint<F>;

        fn domain(&self) -> MontgomeryCurve<F> {
            self.domain
        }

        fn codomain(&self) -> Self::CodomainCurve {
            self.codomain
        }

        fn degree(&self) -> u64 {
            Self::degree(&self.domain, &self.pt_deg.pt)
        }

        fn evaluate(&self, p: &Self::DomainPoint) -> Self::CodomainPoint {
            let x_pts = kernel_rge_prime_degree(&self.domain, &self.pt_deg);
            evaluate(&x_pts, p)
        }

        fn is_separable(&self) -> Choice {
            todo!()
        }
    }
}