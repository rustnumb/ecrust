<h1 align="center">ecrust</h1>
<div align="center">
    <a href="https://github.com/rustnumb/ecrust/actions/workflows/rust.yaml"><img src="https://img.shields.io/github/check-runs/rustnumb/ecrust/main?label=tests&logo=github" alt="GitHub branch checks state"></a>
    |
    <a href="https://crates.io/crates/ecrust"><img src="https://img.shields.io/crates/v/ecrust.svg" alt="Crates.io"></a>
</div>

<p align="center">
  <img src="assets/ecrust_mascot_small.png" alt="ecrust mascot logo" width="420" />
</p>

A Rust workspace for finite fields, elliptic curves, isogenies, and higher-level elliptic-curve protocols.

## Published crates

The project is organized as a layered set of crates:

- `ecrust-fp` — finite-field arithmetic
- `ecrust-ec` — elliptic-curve abstractions and point arithmetic
- `ecrust-isogeny` — isogeny abstractions and scaffolding
- `ecrust-protocol` — protocol-oriented building blocks
- `ecrust` — umbrella crate re-exporting the full stack

Conceptually, the stack looks like this:

```text
protocol       ← ECDH / EC-ElGamal helpers built on curve points
  ↓
isogeny        ← isogeny and kernel abstractions (work in progress)
  ↓
ec             ← elliptic-curve models and point arithmetic
  ↓
fp             ← finite fields: F_p, F_{p^m}, F_2, F_{2^m}
  ↓
crypto-bigint  ← multi-precision integers / Montgomery arithmetic
```

## Recommended dependency choice

If you are starting fresh, the simplest option is the umbrella crate:

```toml
[dependencies]
ecrust = "0.1"
```

This lets you import modules as:

```rust
use ecrust::fp;
use ecrust::ec;
```

If you want finer-grained dependencies, use the layer-specific crates directly. Each crate README shows the recommended `Cargo.toml` entry.

## Current status

The workspace is usable for experiments, API exploration, and teaching, with the following caveats:

- `ecrust-fp` is currently the most complete layer.
- `ecrust-ec` contains several curve models and point types, but parts of the arithmetic are still evolving.
- `ecrust-isogeny` is scaffolding for future work.
- `ecrust-protocol` provides useful building blocks, but should not yet be treated as production-grade cryptographic software.

## Build and test

```bash
cargo build --workspace
cargo test --workspace
```

## Examples and demos

See [DEMO.md](DEMO.md) for concrete examples covering:

- `FieldOps` with `FpElement`
- `FieldOps` with `FpExt`
- `FieldOps` with `F2Ext`
- curve and point instantiation
- `SecretScalar`, `Ecdh`, and `EcElGamal`
- generic helper functions over traits

## Repository layout

```text
ecrust/
├── Cargo.toml
├── README.md
├── DEMO.md
├── ecrust/
│   ├── Cargo.toml
│   ├── README.md
│   └── src/
├── fp/
├── ec/
├── isogeny/
└── protocol/
```

## Disclaimer

This software is currently in **alpha**. We are actively working toward cleaner APIs and broader constant-time coverage, but the code should still be treated as experimental rather than production-hardened.

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache-2.0. See [LICENSE](LICENSE).
