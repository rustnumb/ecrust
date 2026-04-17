# ecrust-ec

Elliptic-curve abstractions and point arithmetic for the **ecrust** ecosystem.

This crate builds on `ecrust-fp` and provides reusable curve and point APIs together with concrete models such as Weierstrass, Montgomery, Edwards, Hessian, and Jacobi forms.

## What is in this crate?

- `CurveOps`: trait for curve-level operations
- `PointOps`: trait for generic group operations on points
- concrete curve and point types for several models
- scalar multiplication helpers and model-specific arithmetic

## Adding it to your project

```toml
[dependencies]
ec = { package = "ecrust-ec", version = "0.1" }
fp = { package = "ecrust-fp", version = "0.1" }
```

Or, with the umbrella crate:

```toml
[dependencies]
ecrust = "0.1"
```

and then import from `ecrust::ec` and `ecrust::fp`.

## Minimal example

```rust
use crypto_bigint::{Uint, const_prime_monty_params};
use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_weierstrass::AffinePoint;
use fp::fp_element::FpElement;

const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

let curve = WeierstrassCurve::new_short(F19::from_u64(2), F19::from_u64(3));
let p = AffinePoint::new(F19::from_u64(1), F19::from_u64(5));

assert!(curve.contains(&p.x, &p.y));
let q = p.double(&curve);
let _r = p.add(&q, &curve);
```

## Status

This crate is suitable for experimentation, education, and API exploration. It is **not yet a hardened production implementation**. In particular, some backends still contain exceptional-case branching and should not be assumed to provide full side-channel resistance.

## Related crates

- `ecrust-fp`: finite-field arithmetic
- `ecrust-isogeny`: isogeny scaffolding on top of this layer
- `ecrust-protocol`: higher-level protocols
- `ecrust`: umbrella crate re-exporting the full stack

## Authors

- Gustavo Banegas
- Martin Azon
- Sam Frengley

## License

Apache-2.0
