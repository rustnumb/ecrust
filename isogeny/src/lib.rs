//! Isogeny computations between elliptic curves.
//!
//! This crate provides algorithms for computing isogenies, their kernels,
//! and the induced maps on points.  It is built on top of:
//!   * `fp`  – finite field arithmetic
//!   * `ec`  – elliptic curve group law

pub mod isogeny;
pub mod kernel;
