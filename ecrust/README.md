# ecrust

The umbrella crate for the **ecrust** ecosystem.

This crate re-exports the workspace layers behind one package so users can start with a single dependency and grow into the modular crates later if they want finer control.

## Re-exported crates

- `ecrust::fp`
- `ecrust::ec`
- `ecrust::isogeny`
- `ecrust::protocol`

## Adding it to your project

```toml
[dependencies]
ecrust = "0.1"
```

## Minimal example

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use ecrust::fp::field_ops::FieldOps;
use ecrust::fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

let a = F19::from_u64(5);
let b = F19::from_u64(9);
assert_eq!((a + b).as_limbs()[0], 14);
assert!(a.invert().into_option().is_some());
```

## Feature flags

All layers are enabled by default.

- `fp`
- `ec`
- `isogeny`
- `protocol`

You can disable default features and select only what you need:

```toml
[dependencies]
ecrust = { version = "0.1", default-features = false, features = ["fp", "ec"] }
```

## When should I use this crate?

Use `ecrust` if you want:

- one dependency for the whole stack
- a simpler onboarding story
- the ability to access modules as `ecrust::fp`, `ecrust::ec`, and so on

Use the individual crates if you want tighter dependencies or if you are building only on one layer.

## Status

The workspace is currently in **alpha**. APIs may change, and the implementation should be considered experimental.

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache-2.0
