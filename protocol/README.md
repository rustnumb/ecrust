# ecrust-protocol

Protocol-oriented building blocks built on top of the **ecrust** field and curve crates.

This crate currently contains small, reusable components such as fixed-width secret scalars, elliptic-curve Diffie-Hellman, and elliptic-curve ElGamal.

## What is in this crate?

- `SecretScalar<LIMBS>`
- `Ecdh`
- `EcElGamal`

## Adding it to your project

```toml
[dependencies]
protocol = { package = "ecrust-protocol", version = "0.1" }
ec = { package = "ecrust-ec", version = "0.1" }
fp = { package = "ecrust-fp", version = "0.1" }
```

Or use the umbrella crate:

```toml
[dependencies]
ecrust = "0.1"
```

and import from `ecrust::protocol`.

## Status

The protocol layer provides a convenient structure for experimentation and teaching, but it should still be treated as **alpha**. It is not yet a production-grade cryptographic toolkit.

## Side-channel note

The crate is designed with constant-time-friendly APIs in mind, but it inherits the security properties of the underlying curve implementations. Until the curve backends are fully hardened, this crate should be considered experimental.

## Related crates

- `ecrust-fp`: finite-field arithmetic
- `ecrust-ec`: elliptic-curve abstractions
- `ecrust-isogeny`: isogeny scaffolding
- `ecrust`: umbrella crate re-exporting the full stack

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache-2.0
