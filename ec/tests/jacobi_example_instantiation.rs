//! Example wrapper showing how to instantiate the reusable suites.
//!
//! Copy this file to your crate's `tests/` directory, then replace the adapter.

#[path = "jacobi_intersection_vectors.rs"]
mod jacobi_intersection_vectors;
#[path = "jacobi_quartic_vectors.rs"]
mod jacobi_quartic_vectors;

// Replace this line with your real adapter file.
// #[path = "field_adapter.rs"]
// mod field_adapter;

// use field_adapter::{f, F};

// jacobi_quartic_test_suite!(ec, F, f);
// jacobi_intersection_test_suite!(ec, F, f);
