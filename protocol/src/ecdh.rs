//! Elliptic-curve Diffie–Hellman helpers.

use ec::point_ops::{PointOps};

use crate::scalar::SecretScalar;

/// Deterministic key pair built from a caller-supplied secret scalar.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KeyPair<P, const LIMBS: usize> {
    pub secret: SecretScalar<LIMBS>,
    pub public: P,
}

/// Stateless ECDH helper.
pub struct Ecdh;

impl Ecdh {
    /// Derive a public key `[secret]G` from a fixed generator/base point.
    pub fn derive_public_key<P, const LIMBS: usize>(
        base_point: &P,
        secret: &SecretScalar<LIMBS>,
        curve: &P::Curve,
    ) -> P
    where
        P: PointOps,
    {
        base_point.scalar_mul(secret.as_limbs(), curve)
    }

    /// Build a deterministic key pair from a secret scalar.
    pub fn keypair_from_secret<P, const LIMBS: usize>(
        base_point: &P,
        secret: SecretScalar<LIMBS>,
        curve: &P::Curve,
    ) -> KeyPair<P, LIMBS>
    where
        P: PointOps,
    {
        let public = Self::derive_public_key(base_point, &secret, curve);
        KeyPair { secret, public }
    }

    /// Compute a shared point `[my_secret] peer_public`.
    ///
    /// Returns `None` when the peer public key is the identity.
    pub fn shared_secret<P, const LIMBS: usize>(
        my_secret: &SecretScalar<LIMBS>,
        peer_public: &P,
        curve: &P::Curve,
    ) -> Option<P>
    where
        P: PointOps,
    {
        if peer_public.is_identity() {
            return None;
        }
        Some(peer_public.scalar_mul(my_secret.as_limbs(), curve))
    }
}
