//! Finite field arithmetic for the `ecrust` library.
//!
//! # Module layout
//!
//! ```text
//! fp
//! ├── f2_element.rs
//! ├── f2_ext.rs
//! ├── field_ops.rs
//! ├── fp_element.rs
//! ├── fp_ext.rs
//! └── lib.rs
//! ```

/// Binary base field $\mathbb{F}_2$ and its arithmetic.
pub mod f2_element;
/// Binary extension fields $\mathbb{F}\_{2^m}$ built from irreducible polynomials.
pub mod f2_ext;
/// Core field traits shared by prime and binary fields.
pub mod field_ops;
/// Prime-field elements over `crypto-bigint` Montgomery arithmetic.
pub mod fp_element;
/// Prime-field extension towers and related helper traits.
pub mod fp_ext;
