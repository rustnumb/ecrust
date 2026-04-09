//! Elliptic curve group operations.
//!
//! This crate provides elliptic curve point arithmetic built on top of
//! the finite field primitives in `fp`.
//!
//! # Module layout
//!
//! ```text
//! ec
//! ├── curve_ops          – Generic Curve trait
//! ├── point_ops          – Generic PointOps trait
//! ├── curve_weierstrass  – General/short Weierstrass curves
//! ├── point_weierstrass  – Affine points on Weierstrass curves
//! ├── curve_montgomery   – Montgomery curves  (By² = x³ + Ax² + x)
//! ├── point_montgomery   – Kummer-line x-only points
//! ├── curve_edwards      – Edwards curves (odd & binary char)
//! └── point_edwards      – Affine points on Edwards curves
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
pub mod curve_montgomery;
pub mod point_montgomery;
pub mod curve_edwards;
pub mod point_edwards;
