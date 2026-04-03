use core::ops::{Add, Mul, Neg, Sub};
use subtle::{ConditionallySelectable, ConstantTimeEq};
use fp::field_ops::FieldOps;
use fp::f2::F2Element;

#[cfg(test)]
mod tests {
    use super::*;

    #[test] fn zero_is_zero() { assert!(F2Element::ZERO.is_zero()); }
    #[test] fn one_is_one() { assert!(F2Element::ONE.is_one()); }
    #[test] fn degree_is_one() { assert_eq!(F2Element::degree(), 1); }
    #[test] fn characteristic_is_2() { assert_eq!(F2Element::characteristic(), vec![2u64]); }
    #[test] fn add_mod_2() {
        // 0 + 1 ≡ 1 (mod 2)
        // 1 + 1 ≡ 0 (mod 2)
        assert_eq!(F2Element::from_u64(0) + F2Element::from_u64(1), F2Element::ONE);
        assert_eq!(F2Element::from_u64(1) + F2Element::from_u64(1), F2Element::ZERO);
    }
    #[test] fn mul_mod_2() {
        // 0 * 0 ≡ 0 (mod 2)
        // 0 * 1 ≡ 0 (mod 2)
        // 1 * 1 ≡ 1 (mod 2)
        assert_eq!(F2Element::from_u64(0) * F2Element::from_u64(0), F2Element::ZERO);
        assert_eq!(F2Element::from_u64(0) * F2Element::from_u64(1), F2Element::ZERO);
        assert_eq!(F2Element::from_u64(1) * F2Element::from_u64(1), F2Element::ONE);
    }
    #[test] fn neg_mod_2()   { assert_eq!(F2Element::ONE.negate(), F2Element::ONE);}
    #[test] fn inv_mod_2()   { assert_eq!(F2Element::ONE.invert().unwrap(), F2Element::ONE);}

}