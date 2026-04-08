# DEMO

This file collects small examples showing how to use the main abstractions in `ecrust`.

A useful mental model is:
- you do **not** instantiate a trait directly;
- you instantiate a **concrete type** that implements the trait.

For example:
- `FpElement<...>` implements `FieldOps`
- `FpExt<...>` implements `FieldOps`
- `F2Ext<...>` implements `FieldOps`
- `WeierstrassCurve<F>` implements `Curve`
- `AffinePoint<F>` implements `PointOps`

---

## 1. Prime field example: instantiate `FieldOps` with `FpElement`

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn main() {
    let a = F19::from_u64(7);
    let b = F19::from_u64(8);

    let add = a + b;                 // 7 + 8 = 15 mod 19
    let mul = a * b;                 // 56 mod 19 = 18
    let sqr = a.square();            // 49 mod 19 = 11
    let inv = a.invert().into_option().unwrap();

    assert_eq!(add.as_limbs()[0], 15);
    assert_eq!(mul.as_limbs()[0], 18);
    assert_eq!(sqr.as_limbs()[0], 11);
    assert!((a * inv).is_one().into());
}
```

Use this pattern whenever you want a concrete prime field `Fp`.

---

## 2. Extension field example: instantiate `FieldOps` with `FpExt`

To define `Fp^m`, you provide an irreducible polynomial by implementing `IrreduciblePoly`.

Here we build `F19^2 = F19[u]/(u^2 + 1)`.

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;
use fp::fp_ext::{FpExt, IrreduciblePoly};

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

struct QuadPoly;
impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
    fn modulus() -> [F19; 2] {
        [F19::one(), F19::zero()] // u^2 + 1
    }
}

type F19_2 = FpExt<Fp19Mod, 1, 2, QuadPoly>;

fn main() {
    let a = F19_2::new([F19::from_u64(3), F19::from_u64(2)]); // 3 + 2u
    let b = F19_2::new([F19::from_u64(5), F19::from_u64(1)]); // 5 + u

    let c = a + b;
    let d = a * b;
    let frob = a.frobenius();
    let norm = a.norm();
    let trace = a.trace();

    let inv = a.invert().into_option().unwrap();
    assert!(FieldOps::mul(&a, &inv).is_one().into());

    let _ = (c, d, frob, norm, trace);
}
```

Use this pattern whenever you want a generic extension field over an odd prime field.

---

## 3. Binary field example: instantiate `FieldOps` with `F2Ext`

To define `F2^m`, implement `BinaryIrreducible` with the full polynomial bitmask.

Here we build `F2^4 = F2[t]/(t^4 + t + 1)`.

```rust
use crypto_bigint::Uint;
use fp::f2_ext::{BinaryIrreducible, F2Ext};
use fp::field_ops::FieldOps;

struct Poly2_4;
impl BinaryIrreducible<1> for Poly2_4 {
    fn modulus() -> Uint<1> { Uint::from_u64(0x13) } // x^4 + x + 1
    fn degree() -> usize { 4 }
}

type GF16 = F2Ext<1, Poly2_4>;

fn main() {
    let a = GF16::from_u64(0b0011);
    let b = GF16::from_u64(0b0101);

    let add = a + b;
    let mul = a * b;
    let sqr = a.square();

    let _ = (add, mul, sqr);
}
```

Use this when you want characteristic-2 arithmetic.

---

## 4. Curve example: instantiate `Curve` with `WeierstrassCurve<F>`

`WeierstrassCurve<F>` is the current concrete curve model in the project.
It implements the generic `Curve` trait.

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use ec::curve_ops::Curve;
use ec::curve_weierstrass::WeierstrassCurve;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(x: u64) -> F19 { F19::from_u64(x) }

fn main() {
    let curve = WeierstrassCurve::new_short(fp(2), fp(3));

    let ainvs = curve.a_invariants();
    let identity = curve.identity();

    let _ = (ainvs, identity);
}
```

For binary fields, you can also use the general constructor `WeierstrassCurve::new(a1, a2, a3, a4, a6)`.

---

## 5. Point example: instantiate `PointOps` with `AffinePoint<F>`

`AffinePoint<F>` is the current concrete point representation and implements `PointOps`.

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_ops::PointOps;
use ec::point_weierstrass::AffinePoint;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(x: u64) -> F19 { F19::from_u64(x) }

fn main() {
    let curve = WeierstrassCurve::new_short(fp(2), fp(3));

    // (1,5) is on y^2 = x^3 + 2x + 3 over F19 because 5^2 = 25 ≡ 6
    // and 1^3 + 2*1 + 3 = 6.
    let p = AffinePoint::new(fp(1), fp(5));
    assert!(curve.contains(&p.x, &p.y));

    let q = p.double(&curve);
    let r = p.add(&q, &curve);
    let s = p.scalar_mul(&[5], &curve);
    let neg = p.negate(&curve);

    assert!(p.add(&neg, &curve).is_identity());
    let _ = (q, r, s);
}
```

---

## 6. Protocol example: instantiate `SecretScalar`, then use `Ecdh`

The protocol layer is generic over the point type. In practice, you choose:
- a concrete field,
- a concrete curve,
- a concrete point type,
- and a fixed-width `SecretScalar`.

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_weierstrass::AffinePoint;
use fp::fp_element::FpElement;
use protocol::ecdh::Ecdh;
use protocol::scalar::SecretScalar;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(x: u64) -> F19 { F19::from_u64(x) }

fn main() {
    let curve = WeierstrassCurve::new_short(fp(2), fp(3));
    let base = AffinePoint::new(fp(1), fp(5));

    let alice_sk = SecretScalar::<1>::new([5]);
    let bob_sk = SecretScalar::<1>::new([7]);

    let alice_pk = Ecdh::derive_public_key(&base, &alice_sk, &curve);
    let bob_pk = Ecdh::derive_public_key(&base, &bob_sk, &curve);

    let alice_ss = Ecdh::shared_secret(&alice_sk, &bob_pk, &curve).unwrap();
    let bob_ss = Ecdh::shared_secret(&bob_sk, &alice_pk, &curve).unwrap();

    assert_eq!(alice_ss, bob_ss);
}
```

This is a toy example to demonstrate the API shape, not a secure deployment choice.

---

## 7. Protocol example: EC-ElGamal on curve points

`EcElGamal` currently encrypts a point message, not arbitrary bytes.

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_weierstrass::AffinePoint;
use fp::fp_element::FpElement;
use protocol::ecdh::Ecdh;
use protocol::elgamal::EcElGamal;
use protocol::scalar::SecretScalar;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(x: u64) -> F19 { F19::from_u64(x) }

fn main() {
    let curve = WeierstrassCurve::new_short(fp(2), fp(3));
    let base = AffinePoint::new(fp(1), fp(5));

    let recipient_sk = SecretScalar::<1>::new([9]);
    let recipient_pk = Ecdh::derive_public_key(&base, &recipient_sk, &curve);

    let message = base.scalar_mul(&[4], &curve);
    let nonce = SecretScalar::<1>::new([3]);

    let ciphertext = EcElGamal::encrypt(&base, &recipient_pk, &message, &nonce, &curve)
        .unwrap();
    let decrypted = EcElGamal::decrypt(&recipient_sk, &ciphertext, &curve);

    assert_eq!(decrypted, message);
}
```

---

## 8. Generic helper over `FieldOps`

One of the main advantages of the trait design is writing code once for several field types.

```rust
use fp::field_ops::FieldOps;

fn eval_monic_quadratic<F: FieldOps>(x: &F, a: &F, b: &F) -> F {
    let x2 = <F as FieldOps>::square(x);
    let ax = <F as FieldOps>::mul(a, x);
    let tmp = <F as FieldOps>::add(&x2, &ax);
    <F as FieldOps>::add(&tmp, b)
}
```

This same helper can be used with:
- `FpElement<...>`
- `FpExt<...>`
- `F2Element`
- `F2Ext<...>`

---

## 9. Generic helper over `PointOps`

Likewise, you can write point code against the trait instead of a concrete point type.

```rust
use ec::point_ops::PointOps;

fn linear_combination<P: PointOps>(p: &P, q: &P, curve: &P::Curve) -> P {
    let two_p = p.double(curve);
    two_p.add(q, curve)
}
```

As long as your point type implements `PointOps`, this helper works unchanged.

---

## 10. Generic helper over `Curve`

If you need the identity point or an on-curve check through the abstract interface:

```rust
use ec::curve_ops::Curve;

fn curve_identity<C: Curve>(curve: &C) -> C::Point {
    curve.identity()
}
```

---

## Notes

- The examples above use very small toy parameters where convenient, because they make the API easier to read.
- For realistic protocol examples, see `protocol/tests/ec_protocols_tests.rs`.
- The current `ec` layer uses affine formulas and still has some exceptional-case branching, so treat these examples as developer documentation and API demonstrations.
- The `isogeny` crate currently exposes the main structs, but evaluation formulas are still TODO.
