# ecRust

<p align="center">
  <img src="assets/ecrust_mascot.png" alt="ecRust mascot logo" width="420" />
</p>

A Rust library for finite-field arithmetic, elliptic-curve operations, isogeny scaffolding, and higher-level elliptic-curve protocols.

The project is organized as a layered set of crates:

```text
protocol       ← ECDH / EC-ElGamal helpers built on curve points
  ↓
isogeny        ← isogeny and kernel abstractions (work in progress)
  ↓
ec             ← elliptic-curve models and point arithmetic
  ↓
fp             ← finite fields: Fp, Fp^m, F2, F2^m
  ↓
crypto-bigint  ← multi-precision integers / Montgomery arithmetic
```

## Workspace crates

### `fp`
Finite-field arithmetic.

Main building blocks:
- `FieldOps`: common trait implemented by field elements.
- `FpElement<MOD, LIMBS>`: prime-field elements over `Fp`.
- `FpExt<MOD, LIMBS, M, P>`: extension-field elements over `Fp^M`.
- `F2Element`: the prime field `F2`.
- `F2Ext<LIMBS, P>`: binary extension fields `F2^m`.
- `IrreduciblePoly` / `BinaryIrreducible`: marker traits used to define extension fields.

### `ec`
Elliptic-curve abstractions and affine Weierstrass arithmetic.

Main building blocks:
- `CurveOps`: generic curve-model trait.
- `PointOps`: generic point/group API.

### `isogeny`
Kernel and isogeny structs.

Current status:
- `KernelSubgroup<C>` exists.
- `Isogeny<C>` exists as the main abstraction.
- evaluation formulas are still TODO.

### `protocol`
Small protocol layer on top of `ec`.

Current modules:
- `SecretScalar<LIMBS>`
- `Ecdh`
- `EcElGamal`

## Current status

This workspace is usable for experiments and API exploration, with the following caveats:
- `fp` is the most complete and best-tested layer.
- `ec` supports affine Weierstrass arithmetic and scalar multiplication, but some methods still contain exceptional-case branching and should not yet be treated as hardened production code.
- `isogeny` is currently scaffolding.
- protocol examples are functional API examples, not production-ready constructions.

## Build and test

```bash
cargo build --workspace
cargo test --workspace
```

## Quick start

### 1. Instantiate a prime field (`FieldOps` via `FpElement`)

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

let a = F19::from_u64(7);
let b = F19::from_u64(8);
let c = a * b;
assert_eq!(c.as_limbs()[0], 18); // 56 mod 19 = 18

let inv = a.invert().into_option().unwrap();
assert!((a * inv).is_one().into());
```

### 2. Instantiate an extension field (`FieldOps` via `FpExt`)

```rust
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;
use fp::fp_ext::{FpExt, IrreduciblePoly};

struct QuadPoly;
impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
    fn modulus() -> [FpElement<Fp19Mod, 1>; 2] {
        [FpElement::one(), FpElement::zero()] // x^2 + 1
    }
}

type F19_2 = FpExt<Fp19Mod, 1, 2, QuadPoly>;

let x = F19_2::new([F19::from_u64(3), F19::from_u64(2)]); // 3 + 2u
let y = x.invert().into_option().unwrap();
assert!(FieldOps::mul(&x, &y).is_one().into());
```

### 3. Instantiate a curve (`Curve`) and a point (`PointOps`)

```rust
use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_weierstrass::AffinePoint;

let curve = WeierstrassCurve::new_short(F19::from_u64(2), F19::from_u64(3));
let p = AffinePoint::new(F19::from_u64(1), F19::from_u64(5));

assert!(curve.contains(&p.x, &p.y));
let q = p.double(&curve);
let r = p.add(&q, &curve);
let s = p.scalar_mul(&[5], &curve);
```

## Examples and demos

See [DEMO.md](DEMO.md) for several concrete examples showing how to instantiate the main traits and concrete types in this workspace:
- `FieldOps` with `FpElement`
- `FieldOps` with `FpExt`
- `FieldOps` with `F2Ext`
- `Curve` / `PointOps` with `WeierstrassCurve` and `AffinePoint`
- `SecretScalar`, `Ecdh`, and `EcElGamal`
- generic helper functions written against traits instead of concrete types

## Repository layout

```text
ecrust/
├── Cargo.toml
├── README.md
├── DEMO.md
├── fp/
│   ├── src/
│   │   ├── field_ops.rs
│   │   ├── fp_element.rs
│   │   ├── fp_ext.rs
│   │   ├── f2_element.rs
│   │   └── f2_ext.rs
│   └── tests/
├── ec/
│   ├── src/
│   │   ├── curve_ops.rs
│   │   ├── point_ops.rs
│   │   ├── curve_weierstrass.rs
│   │   └── point_weierstrass.rs
│   └── tests/
├── isogeny/
│   ├── src/
│   │   ├── kernel.rs
│   │   └── isogeny.rs
│   └── tests/
└── protocol/
    ├── src/
    │   ├── scalar.rs
    │   ├── ecdh.rs
    │   └── elgamal.rs
    └── tests/
```

## Disclaimer
Disclaimer. This software is *currently in an alpha stage*. We are actively working toward constant-time implementations across the project, but achieving this systematically remains an ongoing effort. At this stage, the code should be *treated as experimental*, and it must not be assumed to provide full side-channel resistance or production-grade security guarantees.

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache License 2.0. See [LICENSE](LICENSE).
