//! Elliptic curve group operations.
//!
//! This crate provides elliptic curve point arithmetic built on top of
//! the finite field primitives in `fp`.

extern crate core;

pub mod curve_ops;
pub mod point_ops;

pub mod curve_weierstrass;
pub mod point_weierstrass;
pub mod point_montgomery;
pub mod curve_montgomery;
pub mod point_edwards;
pub mod curve_edwards;
pub mod curve_jacobi_quartic;
pub mod point_jacobi_quartic;
pub mod curve_jacobi_intersection;
pub mod point_jacobi_intersection;
pub mod curve_hessian;
pub mod point_hessian;
/// Point used in the Legendre form.
pub mod point_legendre;
///! Elliptic curve definition in Legendre form.
pub mod curve_legendre;
/// Twisted Hessian curves.
pub mod curve_twisted_hessian;
/// Projective points on twisted Hessian curves.
pub mod point_twisted_hessian;
