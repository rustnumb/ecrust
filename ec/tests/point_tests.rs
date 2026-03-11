//! Integration tests for the `ec` crate.

use ec::point::AffinePoint;

#[test]
fn identity_is_at_infinity() {
    let id = AffinePoint::identity();
    assert!(id.infinity);
}
