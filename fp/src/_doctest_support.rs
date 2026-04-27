//! Hidded modules only used for testing in the documentation

pub mod _doctest_fp_ext {
    use crate::field_ops::FieldOps;
    use crate::fp_element::FpElement;
    use crate::fp_ext::{FpExt, IrreduciblePoly, TonelliShanksConstants};
    use crypto_bigint::{const_prime_monty_params, Uint};

    const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
    pub type Fp19 = FpElement<Fp19Mod, 1>;

    pub struct QuadPoly;
    pub struct TSQuad;

    impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
        fn modulus() -> [Fp19; 2] {
            [Fp19::one(), Fp19::zero()]
        }
    }

    impl TonelliShanksConstants<Fp19Mod, 1, 2, 1> for TSQuad {
        const ORDER: Uint<1> = Uint::<1>::from_u64(360);
        const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(180);
        const S: u64 = 3;
        const T: Uint<1> = Uint::<1>::from_u64(45);
        const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(22);
        const TWOSM1: Uint<1> = Uint::<1>::from_u64(4);
        fn root_of_unity() -> [FpElement<Fp19Mod, 1>; 2] {
            [Fp19::from_u64(3), Fp19::from_u64(3)]
        }
    }

    pub type F19_2 = FpExt<Fp19Mod, 1, 2, 1, QuadPoly, TSQuad>;
}
