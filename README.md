# ecrust

A Rust library for elliptic curve arithmetic and isogeny computations over finite fields, built on top of [`crypto-bigint`](https://crates.io/crates/crypto-bigint).

> **Status:** Early development (v0.1.0). The finite field layer (`fp`) is functional and tested; the elliptic curve (`ec`) and isogeny (`isogeny`) crates contain scaffolding and are under active development.

## Overview

ecrust is organized as a Cargo workspace with three crates that form a layered architecture:

```
isogeny        ← isogeny maps (Vélu's formulas, kernel subgroups)
  ↓
ec             ← elliptic curve point arithmetic (short Weierstrass)
  ↓
fp             ← finite field arithmetic (Fp, Fp^M extensions)
  ↓
crypto-bigint  ← multi-precision integers in Montgomery form
```

### `fp` — Finite Field Arithmetic

The core of the library. Provides constant-time, Montgomery-form field arithmetic for arbitrary primes and their algebraic extensions.

- **`FpElement<MOD, LIMBS>`** — Elements of the prime field Fp = Z/pZ stored in Montgomery representation. Supports addition, subtraction, multiplication, squaring, inversion, square roots, Legendre symbol, and exponentiation via square-and-multiply.
- **`FpExt<MOD, LIMBS, M, P>`** — Elements of the extension field Fp^M = Fp\[x\]/(f(x)) for any degree M and any user-supplied irreducible polynomial. Implements schoolbook polynomial multiplication with modular reduction, inversion via polynomial extended GCD, Frobenius endomorphism, field norm, and field trace.
- **`FieldOps` trait** — A unified algebraic interface that both `FpElement` and `FpExt` implement, enabling generic code over any level of the field tower.
- **`F2Element`** - Elements of the prime field F2 with two elements, stored as Uint<1>. Supports addition, subtraction, multiplication, squaring, inversion, square roots, Legendre symbol, and exponentiation.

### `ec` — Elliptic Curve Group Operations

Defines short Weierstrass curves (y² = x³ + ax + b) and affine point representation. Currently contains type definitions and the point-at-infinity constructor; the full group law is planned.

### `isogeny` — Isogeny Computations

Provides structures for isogeny maps between elliptic curves and their kernel subgroups. Vélu's formulas for isogeny evaluation are planned.

## Getting Started

### Prerequisites

- **Rust** — install via [rustup](https://rustup.rs/)

### Build

```bash
git clone <repo-url>
cd ecrust
cargo build
```

### Test

```bash
cargo test --workspace
```

The `fp` crate includes comprehensive tests covering modular arithmetic over F₁₉, quadratic extension F₁₉² (using x² + 1), and cubic extension F₁₉³ (using x³ − 2).

## Usage

### Defining a prime field

```rust
use crypto_bigint::{const_prime_monty_params, Uint};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

// Define Fp with p = 19
const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

let a = F19::from_u64(7);
let b = F19::from_u64(8);
assert_eq!((a * b).as_limbs()[0], 18); // 56 mod 19 = 18

let inv = a.invert().unwrap();
assert!((a * inv).is_one());
```

### Defining an extension field

```rust
use fp::fp_ext::{FpExt, IrreduciblePoly};

// F₁₉² = F₁₉[x] / (x² + 1)
struct QuadPoly;
impl IrreduciblePoly<Fp19Mod, 1, 2> for QuadPoly {
    fn modulus() -> [FpElement<Fp19Mod, 1>; 2] {
        [FpElement::one(), FpElement::zero()] // x² + 1
    }
}
type F19_2 = FpExt<Fp19Mod, 1, 2, QuadPoly>;

let a = F19_2::new([F19::from_u64(3), F19::from_u64(2)]); // 3 + 2x
let inv = a.invert().unwrap();
assert!(FieldOps::mul(&a, &inv).is_one());

// Frobenius, norm, and trace
let frob = a.frobenius();             // a^p
let n = a.norm();                      // product of conjugates (lies in Fp)
let t = a.trace();                     // sum of conjugates (lies in Fp)
```

## Project Structure

```
ecrust/
├── Cargo.toml              # Workspace root
├── LICENSE                  # Apache 2.0
├── fp/
│   ├── src/
│   │   ├── lib.rs
│   │   ├── field_ops.rs     # FieldOps trait
│   │   ├── fp_element.rs    # Prime field Fp
│   │   └── fp_ext.rs        # Extension field Fp^M
│   └── tests/
│       ├── fp_tests.rs      # F₁₉ tests
│       └── fp_ext_tests.rs  # F₁₉², F₁₉³ tests
├── ec/
│   ├── src/
│   │   ├── lib.rs
│   │   ├── curve.rs         # WeierstrassCurve
│   │   └── point.rs         # AffinePoint
│   └── tests/
│       └── point_tests.rs
└── isogeny/
    ├── src/
    │   ├── lib.rs
    │   ├── isogeny.rs       # Isogeny map (Vélu's formulas)
    │   └── kernel.rs        # KernelSubgroup
    └── tests/
        └── isogeny_tests.rs
```

## Authors 
- Gustavo Banegas 
- Martin Azon
- Sam Frengley

## License

This project is licensed under the Apache License 2.0 — see the [LICENSE](LICENSE) file for details.
