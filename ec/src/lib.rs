//! Elliptic curve group operations.
//!
//! This crate provides elliptic curve point arithmetic built on top of
//! the finite field primitives in `fp`.
//!
//! # Module layout
//!
//! ```text
//! ec
//! ├── curve_ops   – General Curve
//! └── point_ops   – Point and group law (add, double, scalar mul)
//! ```
//!
//! # Supported fields
//!
//! The curve and point types are generic over any `F: FieldOps`, so they
//! work with `FpElement` (prime fields), `FpExt` (prime-field extensions),
//! `F2Element` (binary field), and `F2Ext` (binary-field extensions).

pub mod curve_ops;
pub mod point_ops;

pub mod curve_weierstrass;
pub mod point_weierstrass;
mod point_montgomery;
mod curve_montgomery;
