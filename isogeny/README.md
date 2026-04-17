# ecrust-isogeny

Isogeny abstractions and scaffolding for the **ecrust** ecosystem.

This crate sits above `ecrust-fp` and `ecrust-ec` and is intended to host kernel representations, isogeny structures, and evaluation logic.

## What is in this crate?

- `KernelSubgroup<C>`
- `Isogeny<C>`
- supporting abstractions for future isogeny evaluation code

## Adding it to your project

```toml
[dependencies]
isogeny = { package = "ecrust-isogeny", version = "0.1" }
ec = { package = "ecrust-ec", version = "0.1" }
fp = { package = "ecrust-fp", version = "0.1" }
```

Or depend on the umbrella crate and use `ecrust::isogeny`.

## Status

This crate is currently **work in progress**. The public abstractions are useful for experimentation and for organizing code, but the computational layer is not complete yet.

## Related crates

- `ecrust-fp`: finite-field arithmetic
- `ecrust-ec`: elliptic-curve abstractions
- `ecrust-protocol`: protocol layer
- `ecrust`: umbrella crate re-exporting the full stack

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache-2.0
