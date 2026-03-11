//! Elliptic curve definition (short Weierstrass form).

/// Short Weierstrass curve  y² = x³ + ax + b  over Fp.
pub struct WeierstrassCurve {
    // TODO: replace with fp::FieldElement once defined
    pub a: u64,
    pub b: u64,
}
