//! Protocol layer for `ecrust`.
//!
//! This crate is meant to sit above the reusable finite-field (`fp`) and
//! elliptic-curve (`ec`) crates and host protocol implementations.
//!
//! Current modules:
//! - [`scalar`]   fixed-width secret scalars suitable for constant-time APIs
//! - [`ecdh`]     elliptic-curve Diffie–Hellman key agreement
//! - [`elgamal`]  elliptic-curve ElGamal over group elements
//!
//! # Side-channel note
//!
//! The protocol layer uses a fixed-width scalar container together with the
//! Montgomery-ladder API from `ec::point_ops::CtPointOps`, which keeps the
//! scalar-processing path free of secret-dependent branches.
//!
//! That said, the current affine Weierstrass backend in `ec` still contains
//! exceptional-case branching in point addition and doubling. So these
//! protocols are a good constant-time-oriented structure to build on, but they
//! should not yet be treated as production-grade hardened implementations.

pub mod ecdh;
pub mod elgamal;
pub mod scalar;
