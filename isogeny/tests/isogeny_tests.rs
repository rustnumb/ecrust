//! Integration tests for the `isogeny` crate.

use isogeny::kernel::KernelSubgroup;

#[test]
fn trivial_kernel_contains_identity() {
    let k = KernelSubgroup::trivial();
    assert_eq!(k.points.len(), 1);
    assert!(k.points[0].infinity);
}
