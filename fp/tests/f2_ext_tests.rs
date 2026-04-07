use crypto_bigint::Uint;
use fp::field_ops::FieldOps;
use fp::f2_element::F2Element;
use fp::f2_ext::{F2Ext, BinaryIrreducible};  // ← was missing IrreduciblePoly

struct Deg8Poly;

impl BinaryIrreducible<1> for Deg8Poly {
    fn modulus() -> Uint<1> {
        Uint::<1>::from_u64(0b1_0001_1011)      // x^8 + x^4 + x^3 + x + 1
    }

    fn degree() ->  usize { 8usize }
}

type F256 = F2Ext<1, Deg8Poly>;

fn elt(binary_str: u64) -> F256 { F256::from_u64(binary_str)}


// -----------------------------------------------------------------------
// Structural
// -----------------------------------------------------------------------

#[test]
fn degree_is_8() {
    assert_eq!(F256::degree(), 8usize);
}

#[test]
fn characteristic_equals_() {
    assert_eq!(F256::characteristic(), F2Element::characteristic());
    assert_eq!(F256::characteristic(), vec![2u64]);
}


// -----------------------------------------------------------------------
// Identity elements
// -----------------------------------------------------------------------

#[test]
fn zero_is_zero() {
    assert!(bool::from(F256::zero().is_zero()));
}

#[test]
fn one_is_one() {
    assert!(bool::from(F256::one().is_one()));
}

#[test]
fn nonzero_is_not_zero() {
    assert!(bool::from(!elt(0b000000001).is_zero()));
    assert!(bool::from(!elt(0b1000110).is_zero()));
    assert!(bool::from(!elt(0b1101111010111101).is_zero()));
}

#[test]
fn inv_zero_is_none() {
    assert!(bool::from(F256::zero().invert().is_none()));
}



// -----------------------------------------------------------------------
// Addition / subtraction / negation
// -----------------------------------------------------------------------

#[test]
fn add_coefficient_wise() {
    let a = elt(0b1_1111);     // x^4 + x^3 + x^2 + x + 1
    let b = elt(0b101_0111);   // x^6 + x^4 + x^2 + x + 1
    let c = FieldOps::add(&a, &b);      // x^6 + x^3 = 0b100_1000
    assert_eq!(c, elt(0b100_1000));

    let d = elt(0b10_1010);    // x^5 + x^3 + x
    let e = elt(0b1101_0100);  // x^7 + x^6 + x^4 + x^2
    let f = FieldOps::add(&d, &e);      // x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x = 0b1111_1110
    assert_eq!(f, elt(0b1111_1110));
}

#[test]
fn negate_then_add_is_zero() {
    let a = elt(0b1_1111);     // 1 + x + x^2 + x^3 + x^4
    assert!(bool::from(FieldOps::add(&a, &a.negate()).is_zero()));
}

#[test]
fn add_zero_is_identity() {
    let a = elt(0b1_1111);     // x^4 + x^3 + x^2 + x + 1
    assert_eq!(FieldOps::add(&a, &F256::zero()), a);
}

#[test]
fn double_equals_zero() {
    let a = elt(0b1_1111);     // x^4 + x^3 + x^2 + x + 1
    assert!(bool::from(a.double().is_zero()));
}


// -----------------------------------------------------------------------
// Multiplication and reduction
// -----------------------------------------------------------------------

#[test]
fn mul_by_one_is_identity() {
    let a = elt(0b1_1111);     // x^4 + x^3 + x^2 + x + 1
    assert_eq!(FieldOps::mul(&a, &F256::one()), a);
}


#[test]
fn reduce_concrete_values() {
    let a = elt(0b1_0000_0000);     // x^8
    let b = elt(0b100_0010_0000);   // x^10 + x^5
    let c = elt(0b1_1011);          // x^4 + x^3 + x + 1 = 0b1_1011
    let d = elt(0b100_1100);         // x^6 + x^3 + x^2 = 0b100_1100
    assert_eq!(a, c);
    assert_eq!(b, d);
}


#[test]
fn mul_concrete_values() {
    let a = elt(0b1_1001);        // x^4 + x^3 + 1
    let b = elt(0b101_0011);      // x^6 + x^4 + x + 1
    let c = FieldOps::mul(&a, &b);         // x^7 + x^5 + x^4 + x^3 + x = 0b1011_1010
    assert_eq!(c, elt(0b1011_1010));


    let d = elt(0b10_1010);    // x^5 + x^3 + x
    let e = elt(0b1101_0100);  // x^7 + x^6 + x^4 + x^2
    let f = FieldOps::mul(&d, &e);      // x^7 + x^3 + x^2 = 0b1000_1100
    assert_eq!(f, elt(0b1000_1100));

    let g = FieldOps::mul(&b, &d);      // x^7 + x^2 + x = 0b1000_0110
    assert_eq!(g, elt(0b1000_0110));
}

#[test]
fn mul_commutativity() {
    let a = elt(0b1_1001);        // x^4 + x^3 + 1
    let b = elt(0b101_0011);      // x^6 + x^4 + x + 1
    assert_eq!(FieldOps::mul(&a, &b), FieldOps::mul(&b, &a));
}

#[test]
fn mul_associativity() {
    let a = elt(0b1_1001);        // x^4 + x^3 + 1
    let b = elt(0b101_0011);      // x^6 + x^4 + x + 1
    let c = FieldOps::mul(&a, &b);         // x^7 + x^5 + x^4 + x^3 + x = 0b1011_1010
    assert_eq!(
        FieldOps::mul(&FieldOps::mul(&a, &b), &c),
        FieldOps::mul(&a, &FieldOps::mul(&b, &c))
    );
}

#[test]
fn mul_distributivity() {
    let a = elt(0b1_1001);        // x^4 + x^3 + 1
    let b = elt(0b101_0011);      // x^6 + x^4 + x + 1
    let c = FieldOps::mul(&a, &b);         // x^7 + x^5 + x^4 + x^3 + x = 0b1011_1010
    assert_eq!(
        FieldOps::mul(&a, &FieldOps::add(&b, &c)),
        FieldOps::add(&FieldOps::mul(&a, &b), &FieldOps::mul(&a, &c))
    );
}

#[test]
fn square_equals_mul_self() {
    let a = elt(0b1_1001);        // x^4 + x^3 + 1
    assert_eq!(a.square(), FieldOps::mul(&a, &a));
}


// -----------------------------------------------------------------------
// Inversion
// -----------------------------------------------------------------------

#[test]
fn invert_concrete_value() {
    // a = x^4 + x^3 + 1,  inv_a = x^5 + x^4 + x^3 + x^2 + x + 1 = 0b11_1111
    let a = elt(0b1_1001);
    let inv_a = a.invert().expect("(x^4 + x^3 + 1) is invertible in F_{256}");
    assert_eq!(inv_a, elt(0b11_1111));

    // b = x^6 + x^4 + x + 1,  inv_b = x^7 + x^6 + x^3 + x = 0b1100_1010
    let b = elt(0b101_0011);
    let inv_b = b.invert().expect("(x^6 + x^4 + x + 1) is invertible in F_{256}");
    assert_eq!(inv_b, elt(0b1100_1010));

    // c = x^7 + x^5 + x^4 + x^3 + x,  inv_c = x^6 + x^5 + x^4 + x^2 + x = 0b111_0110
    let c = FieldOps::mul(&a, &b);
    let inv_c = c.invert().expect("(x^7 + x^5 + x^4 + x^3 + x) is invertible in F_{256}");
    assert_eq!(inv_c, elt(0b111_0110));

    // d = x^5 + x^3 + x, inv_d = x^7 + x^4 + x^3 = 0b1001_1000
    let d = elt(0b10_1010);
    let inv_d = d.invert().expect("(x^5 + x^3 + x) is invertible in F_{256}");
    assert_eq!(inv_d, elt(0b1001_1000));
}

#[test]
fn invert_zero_is_none() {
    assert!(bool::from(F256::zero().invert().is_none()));
}


// -----------------------------------------------------------------------
// Frobenius
// -----------------------------------------------------------------------

#[test]
fn frobenius_concrete_values() {
    let a = elt(0b1_1001);         // x^4 + x^3 + 1
    let a_sq = a.frobenius();               // x^6 + x^4 + x^3 + x = 0b101_1010
    assert_eq!(a_sq, elt(0b101_1010));
    assert_eq!(a_sq.frobenius(), elt(0b1111_0100));   // x^7 + x^6 + x^5 + x^4 + x^2 = 0b1111_0100
}


// -----------------------------------------------------------------------
// Norm and trace
// -----------------------------------------------------------------------

#[test]
fn norm_concrete_values() {
    let a = elt(0b1_1001);         // x^4 + x^3 + 1
    assert_eq!(a.norm(), elt(0b1));
    let b = elt(0b101_0011);       // x^6 + x^4 + x + 1
    assert_eq!(b.norm(), elt(0b1));
}

#[test]
fn trace_concrete_values() {
    let a = elt(0b1_1001);         // x^4 + x^3 + 1
    assert_eq!(a.trace(), elt(0b0));
    let b = elt(0b101_0011);       // x^6 + x^4 + x + 1
    assert_eq!(b.trace(), elt(0b0));
    let e = elt(0b1101_0100);      // x^7 + x^6 + x^4 + x^2
    assert_eq!(e.trace(), elt(0b1));
}


// -----------------------------------------------------------------------
// pow
// -----------------------------------------------------------------------

#[test]
fn pow_zero_is_one()  {
    let a = elt(0b1_1001);         // x^4 + x^3 + 1
    assert!(bool::from(a.pow(&[0]).is_one()));
}

#[test]
fn pow_one_is_self()  {
    let a = elt(0b1_1001);         // x^4 + x^3 + 1
    assert_eq!(a.pow(&[1]), a);
}

#[test]
fn pow_concrete_values() {
    let a = elt(0b1_1001);                            // x^4 + x^3 + 1
    assert_eq!(a.pow(&[39]), elt(0b110_0010));         // a^39 = x^6 + x^5 + x = 0b110_0010
    let b = elt(0b101_0011);                          // x^6 + x^4 + x + 1
    assert_eq!(b.pow(&[73]), elt(0b100_0011));         // b^73 = x^6 + x + 1 = 0b100_0011
    let e = elt(0b1101_0100);                          // x^7 + x^6 + x^4 + x^2
    assert_eq!(e.pow(&[111]), elt(0b1000));             // e^111 = x^3 = 0b1000
}

#[test]

fn pow_group_order() {
    let a = elt(0b1_1001);                            // x^4 + x^3 + 1
    assert_eq!(a.pow(&[255]), elt(0b1));
}



// -----------------------------------------------------------------------
// sqrt
// -----------------------------------------------------------------------


#[test]
fn sqrt_of_zero_is_zero()  {
    let a = F256::zero();         // 0
    assert_eq!(a.sqrt().expect("zero is a square!"), F256::zero());
}

#[test]
fn sqrt_of_one_is_one()  {
    let a = F256::one();         // 1
    assert_eq!(a.sqrt().expect("one is a square!"), F256::one());
}

#[test]
fn sqrt_concrete_values() {
    // a = x^4 + x^3 + 1,  sqrt_a = x^7 + x^6 + x^5 + x^3 + x = 0b1110_1010
    let a = elt(0b1_1001);
    let sqrt_a = a.sqrt().expect("(x^4 + x^3 + 1) is a square in F_{256}");
    assert_eq!(sqrt_a, elt(0b1110_1010));

    // b = x^6 + x^4 + x + 1,  sqrt_b = x^7 + x^6 + x^5 + x^4 + x^2 + x + 1 = 0b1111_0111
    let b = elt(0b101_0011);
    let sqrt_b = b.sqrt().expect("(x^6 + x^4 + x + 1) is a square in F_{256}");
    assert_eq!(sqrt_b, elt(0b1111_0111));

    // d = x^5 + x^3 + x, sqrt_d = x^7 + x^6 + x^4 = 0b1101_0000
    let d = elt(0b10_1010);
    let sqrt_d = d.sqrt().expect("(x^5 + x^3 + x) is a square in F_{256}");
    assert_eq!(sqrt_d, elt(0b1101_0000));

    // e = x^7 + x^6 + x^4 + x^2, sqrt_e = x^7 + x^4 + x^3 + x^2 + x + 1 = 0b1001_1111
    let e = elt(0b1101_0100);
    let sqrt_e = e.sqrt().expect("(x^7 + x^6 + x^4 + x^2) is a square in F_{256}");
    assert_eq!(sqrt_e, elt(0b1001_1111));
}

#[test]
fn sqrt_of_product() {
    // a = x^4 + x^3 + 1,  sqrt_a = x^5 + x^4 + x^3 + x^2 + x + 1 = 0b11_1111
    let a = elt(0b1_1001);
    let sqrt_a = a.sqrt().expect("(x^4 + x^3 + 1) is a square in F_{256}");

    // b = x^6 + x^4 + x + 1,  sqrt_b = x^7 + x^6 + x^3 + x = 0b1100_1010
    let b = elt(0b101_0011);
    let sqrt_b = b.sqrt().expect("(x^6 + x^4 + x + 1) is a square in F_{256}");

    // c = x^7 + x^5 + x^4 + x^3 + x,  sqrt_c = x^6 + x^5 + x^4 + x^2 + x = 0b111_0110
    let c = FieldOps::mul(&a, &b);
    let sqrt_c = c.sqrt().expect("(x^7 + x^5 + x^4 + x^3 + x) is a square in F_{256}");
    assert_eq!(sqrt_c, FieldOps::mul(&sqrt_a, &sqrt_b));
}