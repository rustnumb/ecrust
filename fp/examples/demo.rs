//! Demo: random sampling + pretty printing for every field type.
//!
//! Run with:  cargo run --bin demo
//!
//! Adjust the concrete type aliases at the top if your project already
//! defines its own moduli / irreducible polynomials.

use crypto_bigint::{const_prime_monty_params, Uint};


// Pull in the library types
use fp::f2_element::F2Element;
use fp::f2_ext::{BinaryIrreducible, F2Ext};
use fp::field_ops::{FieldOps, FieldRandom};
use fp::fp_element::FpElement;
use fp::fp_ext::{FpExt, IrreduciblePoly, TonelliShanksConstants};

// ===================================================================
// 1. Small prime field  F_19  (single-limb)
// ===================================================================

// The hex representation of 19 is 0x13, padded to 16 hex digits for
// a single 64-bit limb.
const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);

// Safety: 19 is prime.

type Fp19 = FpElement<Fp19Mod, 1>;

// ===================================================================
// 2. Extension field  F_{19^2} = F_19[x] / (x^2 + 1)
// ===================================================================

struct Poly19Quad;

impl IrreduciblePoly<Fp19Mod, 1, 2> for Poly19Quad {
    // f(x) = x^2 + 1  →  non-leading coeffs [c0, c1] = [1, 0]
    fn modulus() -> [FpElement<Fp19Mod, 1>; 2] {
        [Fp19::from_u64(1), Fp19::from_u64(0)]
    }
}

// Tonelli-Shanks constants for F_{19^2}
// |F_{19^2}*| = 19^2 - 1 = 360 = 2^3 * 45
struct TS19Quad;

impl TonelliShanksConstants<Fp19Mod, 1, 2, 1> for TS19Quad {
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

type Fp19_2 = FpExt<Fp19Mod, 1, 2, 1, Poly19Quad, TS19Quad>;


// ===================================================================
// 5. Prime field  F_97283
// ===================================================================

const_prime_monty_params!(Fp97283Mod, Uint<1>, "0000000000017C03", 2);

type Fp97283 = FpElement<Fp97283Mod, 1>;

// ===================================================================
// 6. Extension field  F_{97283^2} = F_p[x] / (x^2 + 1)
//    Since 97283 ≡ 3 mod 4, x^2 + 1 is irreducible over F_p.
// ===================================================================

struct Poly97283Quad;

impl IrreduciblePoly<Fp97283Mod, 1, 2> for Poly97283Quad {
    fn modulus() -> [FpElement<Fp97283Mod, 1>; 2] {
        // x^2 + 1  -> [1, 0]
        [Fp97283::from_u64(1), Fp97283::from_u64(0)]
    }
}

struct TS97283Quad;

impl TonelliShanksConstants<Fp97283Mod, 1, 2, 1> for TS97283Quad {
    // p^2 - 1 = 9463982088 = 2^3 * 1182997761
    const ORDER: Uint<1> = Uint::<1>::from_u64(9_463_982_088);
    const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(4_731_991_044);
    const S: u64 = 3;
    const T: Uint<1> = Uint::<1>::from_u64(1_182_997_761);
    const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(591_498_880);
    const TWOSM1: Uint<1> = Uint::<1>::from_u64(4);

    fn root_of_unity() -> [FpElement<Fp97283Mod, 1>; 2] {
        // An element of exact order 8 in F_{p^2}
        [Fp97283::from_u64(96_901), Fp97283::from_u64(382)]
    }
}

type Fp97283_2 = FpExt<Fp97283Mod, 1, 2, 1, Poly97283Quad, TS97283Quad>;

// ===================================================================
// 7. Extension field  F_{97283^3} = F_p[x] / (x^3 + x + 1)
// ===================================================================

struct Poly97283Cubic;

impl IrreduciblePoly<Fp97283Mod, 1, 3> for Poly97283Cubic {
    fn modulus() -> [FpElement<Fp97283Mod, 1>; 3] {
        // x^3 + x + 1  -> [1, 1, 0]
        [
            Fp97283::from_u64(1),
            Fp97283::from_u64(1),
            Fp97283::from_u64(0),
        ]
    }
}

struct TS97283Cubic;

impl TonelliShanksConstants<Fp97283Mod, 1, 3, 1> for TS97283Cubic {
    // p^3 - 1 = 920684569564186 = 2^1 * 460342284782093
    const ORDER: Uint<1> = Uint::<1>::from_u64(920_684_569_564_186);
    const HALF_ORDER: Uint<1> = Uint::<1>::from_u64(460_342_284_782_093);
    const S: u64 = 1;
    const T: Uint<1> = Uint::<1>::from_u64(460_342_284_782_093);
    const PROJENATOR_EXP: Uint<1> = Uint::<1>::from_u64(230_171_142_391_046);
    const TWOSM1: Uint<1> = Uint::<1>::from_u64(1);

    fn root_of_unity() -> [FpElement<Fp97283Mod, 1>; 3] {
        // Same convention as your F19^3 example when S = 1
        [
            Fp97283::from_u64(1),
            Fp97283::from_u64(0),
            Fp97283::from_u64(0),
        ]
    }
}

type Fp97283_3 = FpExt<Fp97283Mod, 1, 3, 1, Poly97283Cubic, TS97283Cubic>;

// ===================================================================
// 8. Extension field  F_{97283^4} = F_p[x] / (x^4 + 2x^2 + 2)
//    This one needs N = 2 because p^4 - 1 does not fit in 64 bits.
// ===================================================================

struct Poly97283Quartic;

impl IrreduciblePoly<Fp97283Mod, 1, 4> for Poly97283Quartic {
    fn modulus() -> [FpElement<Fp97283Mod, 1>; 4] {
        // x^4 + 2x^2 + 2  -> [2, 0, 2, 0]
        [
            Fp97283::from_u64(2),
            Fp97283::from_u64(0),
            Fp97283::from_u64(2),
            Fp97283::from_u64(0),
        ]
    }
}

struct TS97283Quartic;

impl TonelliShanksConstants<Fp97283Mod, 1, 4, 2> for TS97283Quartic {
    // p^4 - 1 = 89566956980912803920 = 2^4 * 5597934811307050245
    const ORDER: Uint<2> =
        Uint::<2>::from_words([0xdafdc0d3fc005050, 0x0000000000000004]);

    const HALF_ORDER: Uint<2> =
        Uint::<2>::from_words([0x6d7ee069fe002828, 0x0000000000000002]);

    const S: u64 = 4;

    const T: Uint<2> =
        Uint::<2>::from_words([0x4dafdc0d3fc00505, 0x0000000000000000]);

    const PROJENATOR_EXP: Uint<2> =
        Uint::<2>::from_words([0x26d7ee069fe00282, 0x0000000000000000]);

    const TWOSM1: Uint<2> =
        Uint::<2>::from_words([0x0000000000000008, 0x0000000000000000]);

    fn root_of_unity() -> [FpElement<Fp97283Mod, 1>; 4] {
        // An element of exact order 16 in F_{p^4}
        [
            Fp97283::from_u64(0),
            Fp97283::from_u64(0),
            Fp97283::from_u64(0),
            Fp97283::from_u64(70_178),
        ]
    }
}

type Fp97283_4 = FpExt<Fp97283Mod, 1, 4, 2, Poly97283Quartic, TS97283Quartic>;

// ===================================================================
// 3. Binary extension field  F_{2^8} = F_2[x] / (x^8 + x^4 + x^3 + x + 1)
//    This is the AES / Rijndael polynomial.
// ===================================================================

struct AesPoly;

impl BinaryIrreducible<1> for AesPoly {
    fn modulus() -> Uint<1> {
        // x^8 + x^4 + x^3 + x + 1  =  0b1_0001_1011  =  0x11B
        Uint::<1>::from_u64(0x11B)
    }
    fn degree() -> usize {
        8
    }
}

struct F2511Poly;

impl BinaryIrreducible<8> for F2511Poly {
    fn modulus() -> Uint<8> {
        let one = Uint::<8>::from_u64(1);
        (one << 511) | (one << 10) | one
    }

    fn degree() -> usize {
        511
    }
}

type GF2_511 = F2Ext<8, F2511Poly>;
type GF256 = F2Ext<1, AesPoly>;

// ===================================================================
// main
// ===================================================================

fn main() {
    let mut rng = rand::rng();
    // --- F_2 -----------------------------------------------------------
    println!("===== F_2 (base binary field) =====");
    for i in 0..5 {
        let a = F2Element::random(&mut rng);
        println!("  sample {i}: {a}");
    }

    // --- F_{2^8} (GF(256)) --------------------------------------------
    println!("\n===== F_{{2^8}} (AES field, polynomial form) =====");
    for i in 0..5 {
        let a = GF256::random(&mut rng);
        println!("  sample {i}: {a}");
    }

    // quick sanity: a * a^{-1} == 1
    let a = GF256::random(&mut rng);
    let a_inv = a.invert().unwrap();
    let prod = a * a_inv;
    println!("  sanity  a  = {a}");
    println!("  sanity 1/a = {a_inv}");
    println!("  a * 1/a    = {prod}  (should be 1)");

    println!("\n===== F_{{2^511}} with modulus x^511 + x^10 + 1 =====");
    for i in 0..5 {
        let a = GF2_511::random(&mut rng);
        println!("  sample {i}: {a}");
    }

    let a = GF2_511::random(&mut rng);
    let a_inv = a.invert().unwrap();
    let prod = a * a_inv;
    println!("  sanity  a  = {a}");
    println!("  sanity 1/a = {a_inv}");
    println!("  a * 1/a    = {prod}  (should be 1)");

    // --- F_19 ----------------------------------------------------------
    println!("\n===== F_19 (prime field) =====");
    for i in 0..5 {
        let a = Fp19::random(&mut rng);
        println!("  sample {i}: {a}");
    }

    // --- F_{19^2} ------------------------------------------------------
    println!("\n===== F_{{19^2}} = F_19[x]/(x^2+1) =====");
    for i in 0..5 {
        let a = Fp19_2::random(&mut rng);
        println!("  sample {i}: {a}");
    }

    // quick sanity: a * a^{-1} == 1
    let b = Fp19_2::random(&mut rng);
    let b_inv = b.invert().unwrap();
    let prod = b * b_inv;
    println!("  sanity  b  = {b}");
    println!("  sanity 1/b = {b_inv}");
    println!("  b * 1/b    = {prod}  (should be 1)");

    println!("\n===== F_97283 (prime field) =====");
    for i in 0..5 {
        let a = Fp97283::random(&mut rng);
        println!("  sample {i}: {a}");
    }

    println!("\n===== F_{{97283^2}} = F_97283[x]/(x^2+1) =====");
    for i in 0..5 {
        let a = Fp97283_2::random(&mut rng);
        println!("  sample {i}: {a}");
    }
    let a2 = Fp97283_2::random(&mut rng);
    let a2_inv = a2.invert().unwrap();
    println!("  sanity  a  = {a2}");
    println!("  sanity 1/a = {a2_inv}");
    println!("  a * 1/a    = {}  (should be 1)", a2 * a2_inv);

    println!("\n===== F_{{97283^3}} = F_97283[x]/(x^3+x+1) =====");
    for i in 0..5 {
        let a = Fp97283_3::random(&mut rng);
        println!("  sample {i}: {a}");
    }
    let a3 = Fp97283_3::random(&mut rng);
    let a3_inv = a3.invert().unwrap();
    println!("  sanity  a  = {a3}");
    println!("  sanity 1/a = {a3_inv}");
    println!("  a * 1/a    = {}  (should be 1)", a3 * a3_inv);

    println!("\n===== F_{{97283^4}} = F_97283[x]/(x^4+2x^2+2) =====");
    for i in 0..5 {
        let a = Fp97283_4::random(&mut rng);
        println!("  sample {i}: {a}");
    }
    let a4 = Fp97283_4::random(&mut rng);
    let a4_inv = a4.invert().unwrap();
    println!("  sanity  a  = {a4}");
    println!("  sanity 1/a = {a4_inv}");
    println!("  a * 1/a    = {}  (should be 1)", a4 * a4_inv);

}