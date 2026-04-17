# ecrust-fp

Finite-field arithmetic for the **ecrust** ecosystem.

This crate is the lowest layer of the project and provides reusable field types and traits for prime fields, binary fields, and extension fields.

## What is in this crate?

- `FieldOps`: shared trait for field arithmetic
- `FpElement<MOD, LIMBS>`: elements of a prime field $\mathbb{F}_p$
- `FpExt<...>`: extension fields over prime fields
- `F2Element`: the binary field $\mathbb{F}_2$
- `F2Ext<...>`: binary extension fields $\mathbb{F}_{2^m}$

## Adding it to your project

Because the published package name is `ecrust-fp` while the library crate name is `fp`, the most convenient dependency declaration is:

```toml
[dependencies]
fp = { package = "ecrust-fp", version = "0.1" }
```

If you prefer a single entry point, you can also depend on the umbrella crate:

```toml
[dependencies]
ecrust = "0.1"
```

and then use `ecrust::fp`.

## Minimal example

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use fp::field_ops::FieldOps;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

let a = F19::from_u64(7);
let b = F19::from_u64(8);
let c = a * b;
assert_eq!(c.as_limbs()[0], 18);

let inv = a.invert().into_option().unwrap();
assert!((a * inv).is_one().into());
```

## Status

`ecrust-fp` is currently the most mature layer in the workspace, but the project as a whole is still in **alpha**. APIs may evolve, and the code should be treated as experimental rather than production-hardened.

## Related crates

- `ecrust-ec`: elliptic-curve abstractions built on top of this crate
- `ecrust-isogeny`: isogeny scaffolding built on top of `fp` and `ec`
- `ecrust-protocol`: higher-level elliptic-curve protocols
- `ecrust`: umbrella crate re-exporting the full stack

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache-2.0
