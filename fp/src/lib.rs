#![doc = include_str!("../../katex-header.html")]
//! Finite field arithmetic for the `ecrust` library.
//!
//! # Module layout
//!
//! ```text
//! fp
//! ├── field_ops          - FieldOps trait (the algebraic contract)
//! ├── fp_element         - Base prime field Fp element
//! └── fp_ext             - Extension prime field Elements
//! ```

/// Binary base field $\mathbb{F}_2$ and its arithmetic.
pub mod f2_element;
/// Binary extension fields $\mathbb{F}_{2^m}$ built from irreducible polynomials.
pub mod f2_ext;
/// Core field traits shared by prime and binary fields.
pub mod field_ops;
/// Prime-field elements over `crypto-bigint` Montgomery arithmetic.
pub mod fp_element;
/// Prime-field extension towers and related helper traits.
pub mod fp_ext;
