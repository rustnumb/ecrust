//! Finite field arithmetic for the `ecrust` library.
//!
//! # Module layout
//!
//! ```text
//! fp
//! ├── field_ops          – FieldOps trait (the algebraic contract)
//! ├── fp_element         – Base prime field Fp element
//! └── fp_extension_arithmetic
//! ```

pub mod field_ops;
pub mod fp_element;
