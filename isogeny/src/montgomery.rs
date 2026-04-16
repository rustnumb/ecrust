//! Isogenies for montgomery curves
//!
//! AT THE MOMENT THIS WORKS ONLY IN ODD CHARACTERISTIC!!!


use subtle::{Choice, CtOption, ConditionallySelectable, ConstantTimeEq};

use fp::field_ops::FieldOps;
use ec::curve_montgomery::MontgomeryCurve;
use ec::point_montgomery::KummerPoint;

#[derive(Debug, Clone)]
pub struct MontgomeryIsogeny<F: FieldOps + Copy> {
    pub domain: MontgomeryCurve<F>,
    pub codomain: MontgomeryCurve<F>,
    // x-coordinate of the generator of the kernel <K>
    pub xK: KummerPoint<F>,
}


// ---------------------------------------------------------------------------
// Manual trait impls
// ---------------------------------------------------------------------------

impl<F: FieldOps> PartialEq for MontgomeryIsogeny<F> {
    fn eq(&self, other: &Self) -> bool {
        self.domain == other.domain && self.codomain == other.codomain && self.xK == other.xK
    }
}

impl<F: FieldOps> Eq for MontgomeryIsogeny<F>
where
    F: FieldOps + ConstantTimeEq,
{ }


// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

impl<F: FieldOps> MontgomeryIsogeny<F> {
    pub fn new(domain: MontgomeryCurve<F>, xK: KummerPoint<F>) -> Self {
        let (a_codomain, b_codomain) = Self::new_codomain_helper(&domain, xK);
        Self{
            domain,
            codomain: MontgomeryCurve::new(a_codomain, b_codomain),
            xK
        }
    }

    fn new_codomain_helper(domain: & MontgomeryCurve<F>, xK: KummerPoint<F>) -> (F, F) {
        todo!()
    }
}




// ---------------------------------------------------------------------------
// Constant-time functionalities
// ---------------------------------------------------------------------------













