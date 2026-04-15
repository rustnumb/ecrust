use ec::curve_ops::Curve;
#[path = "jacobi_quartic_vectors.rs"]
mod jacobi_quartic_vectors;

#[path = "jacobi_intersection_vectors.rs"]
mod jacobi_intersection_vectors;

#[path = "field_adapter.rs"]
mod field_adapter;

use field_adapter::{f, F};

// Because the suites use #[macro_export], call them directly:
jacobi_quartic_test_suite!(ec, F, f);
jacobi_intersection_test_suite!(ec, F, f);


