#![doc = include_str!("../../katex-header.html")]
//! Finite field arithmetic for the `ecrust` library.
//!
//! # Module layout
//!
//! ```text
//! fp
//! ├── field_ops          – FieldOps trait (the algebraic contract)
//! ├── fp_element         – Base prime field Fp element
//! └── fp_ext -Extension prime field Elements
//! ```

pub mod f2_element;
pub mod f2_ext;
pub mod field_ops;
pub mod fp_element;
pub mod fp_ext;
