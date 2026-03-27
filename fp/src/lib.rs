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

pub mod field_ops;
pub mod fp_element;
pub mod fp_ext;
pub mod f2;
mod f2m;
