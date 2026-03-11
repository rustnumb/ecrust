//! Integration tests for the `fp` crate.

use fp::add;

#[test]
fn add_two_field_elements() {
    assert_eq!(add(3, 5), 8);
}
