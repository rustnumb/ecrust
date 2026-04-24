/*

//! Isogenies for montgomery curves
//!
//! AT THE MOMENT THIS WORKS ONLY IN ODD CHARACTERISTIC!!!


use subtle::{ConstantTimeEq};

use fp::field_ops::FieldOps;
use ec::curve_montgomery::MontgomeryCurve;
use ec::point_montgomery::KummerPoint;
use fp::{ref_field_impl, ref_field_trait_impl};

#[derive(Debug, Clone)]
pub struct MontgomeryIsogeny<F: FieldOps + Copy> {
    pub domain: MontgomeryCurve<F>,
    pub codomain: MontgomeryCurve<F>,
    // x-coordinate of the generator of the kernel <K>
    pub xk: KummerPoint<F>,
}


// ---------------------------------------------------------------------------
// Manual trait impls
// ---------------------------------------------------------------------------

ref_field_trait_impl! {
    impl<F> PartialEq for MontgomeryIsogeny<F> {
        fn eq(&self, other: &Self) -> bool {
            self.domain == other.domain && self.codomain == other.codomain && self.xk == other.xk
        }
    }
}


impl<F: FieldOps> Eq for MontgomeryIsogeny<F>
where
    F: FieldOps + ConstantTimeEq,
{ }


// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------


ref_field_impl! {
    impl<F> MontgomeryIsogeny<F> {
        pub fn new(domain: MontgomeryCurve<F>, xk: KummerPoint<F>) -> Self {
            let (a_codomain, b_codomain) = Self::new_codomain_helper(&domain, xk);
            Self{
                domain,
                codomain: MontgomeryCurve::new(a_codomain, b_codomain),
                xk
            }
        }

        fn new_codomain_helper(domain: & MontgomeryCurve<F>, xk: KummerPoint<F>) -> (F, F) {
            todo!()
        }
    }
}




// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------













*/
