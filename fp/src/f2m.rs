// Generic finary fields F_{2^m} = F_2[x] / (f(x))

use std::marker::PhantomData;
use crypto_bigint::{Uint, CtSelect};


// ===========================================================================
// IrreduciblePoly — the only thing callers need to implement for a new field
// ===========================================================================

pub trait BinaryIrreducible<const LIMBS: usize>: 'static {
    // Full polynomial bitmask, including the leading term x^m
    fn modulus() -> Uint<LIMBS>;

    // Degree m of the irreducible polynomial
    fn degree() ->  usize;
}


// ===========================================================================
// F2Ext — element of F_{2^M}
// ===========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct F2Ext<const LIMBS: usize, P>
where
    P: BinaryIrreducible<LIMBS>
{
    pub value: Uint<LIMBS>,
    _phantom: PhantomData<P>
}


// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------


impl<const LIMBS: usize, P> F2Ext<LIMBS, P>
where
    P: BinaryIrreducible<LIMBS>
{
    pub fn new(x: Uint<LIMBS>) -> Self {
        Self { value: reduce::<LIMBS, P>(x), _phantom: PhantomData }
    }
    pub fn from_u64(x: u64) -> Self { Self::new(Uint::from(x)) }
    pub fn from_uint(x: Uint<LIMBS>) -> Self { Self::new(x) }
    pub fn as_uint(&self) -> Uint<LIMBS> { self.value }
    pub fn degree() -> usize { P::degree() }


}


// ---------------------------------------------------------------------------
// Private helper functions
// ---------------------------------------------------------------------------

fn reduce<const LIMBS: usize, P>(mut a: Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS> {
    let modulus = P::modulus();
    let m = P::degree();
    let nbits = Uint::<LIMBS>::BITS as usize;

    debug_assert!(m > 0);
    debug_assert!(m < nbits);

    // Fixed iteration count => no leakage from the current degree of `a`.
    for i in (m..nbits).rev() {
        let bit = a.bit(i as u32);      // Choice: whether x^i is present
        let shifted = modulus << (i - m);
        let reduced = a ^ shifted;

        // Constant-time:
        // - keep `a` if bit == 0
        // - replace by `reduced` if bit == 1
        a = Uint::<LIMBS>::ct_select(&a, &reduced, bit.into());
    }
    a
}

// This implements the add-and-shift with immediate modular reduction algorithm
fn mul<const LIMBS: usize, P>(a: Uint<LIMBS>, b: Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>
{
    let m = P::degree();
    let modulus = P::modulus();

    // lower = modulus without the x^m term
    let lower = modulus ^ (Uint::<LIMBS>::ONE << m);

    let mut res = Uint::<LIMBS>::ZERO;
    let mut cur = a;        // current term appearing in the sum

    for i in 0..m {
        let bit = b.bit(i as u32);

        // if bit(b, i) then res += cur (sum here is just xoring)
        let res_xor = res ^ cur;
        res = Uint::<LIMBS>::ct_select(&res, &res_xor, bit.into());

        // multiply cur by x and reduce modulo P
        let top = cur.bit((m-1) as u32);
        let shifted = cur << 1;
        let reduced = shifted ^ lower;      // cur = shifted + top*lower
        cur = Uint::<LIMBS>::ct_select(&shifted, &reduced, top.into());
    }

    res
}


fn square<const LIMBS: usize, P>(a: Uint<LIMBS>) -> Uint<LIMBS>
where
    P: BinaryIrreducible<LIMBS>
{
    let m = P::degree();
    let modulus = P::modulus();

    // lower = modulus without the x^m term
    let lower = modulus ^ (Uint::<LIMBS>::ONE << m);

    let mut res = Uint::<LIMBS>::ZERO;
    let mut cur = Uint::<LIMBS>::ONE;   // cur = x^(2i) mod P, initially i = 0

    for i in 0..m {
        let bit = a.bit(i as u32);

        // if bit(a, i) = 1 then res += x^(2i) mod P
        let res_xor = res ^ cur;
        res = Uint::<LIMBS>::ct_select(&res, &res_xor, bit.into());

        // multiply cur by x^2 and reduce mod P (we do so in 2 steps)
        let top1 = cur.bit((m-1) as u32);
        let shifted1 = cur << 1;
        let reduced1 = shifted1 ^ lower;
        cur = Uint::<LIMBS>::ct_select(&shifted1, &reduced1, top1.into());

        let top2 = cur.bit((m-1) as u32);
        let shifted2 = cur << 1;
        let reduced2 = shifted2 ^ lower;
        cur = Uint::<LIMBS>::ct_select(&shifted2, &reduced2, top2.into());
    }

    res
}
