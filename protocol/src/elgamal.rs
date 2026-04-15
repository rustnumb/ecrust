//! Elliptic-curve ElGamal over group elements.
//!
//! This module deliberately works on *points as messages* instead of trying to
//! encode arbitrary byte strings into curve points. That keeps the abstraction
//! compact and lets it directly reuse the existing `ec` layer.

use ec::point_ops::PointAdd;

use crate::scalar::SecretScalar;

/// ElGamal ciphertext consisting of a pair $(c_1, c_2)$.
///
/// The ciphertext is formed as:
///
/// $$
/// c_1 = rG, \quad c_2 = M + r \mathsf{PK}
/// $$
///
/// where:
/// - $r$ is a random
/// - $G$ is the generator of the group
/// - $PK$ is the recipient’s public key
/// - $M$ is the plaintext message point
///
/// This struct stores the two ciphertext components.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ciphertext<P> {
    /// The ephemeral public component $c_1 = r G$.
    pub ephemeral_public: P,

    /// The blinded message component $c_2 = M + r \mathsf{PK}$.
    pub blinded_message: P,
}

/// Stateless EC-ElGamal helper.
pub struct EcElGamal;

impl EcElGamal {
    /// Encrypt a point `message` using a recipient public key and an
    /// application-provided ephemeral scalar.
    pub fn encrypt<P, const LIMBS: usize>(
        base_point: &P,
        recipient_public: &P,
        message: &P,
        ephemeral_secret: &SecretScalar<LIMBS>,
        curve: &P::Curve,
    ) -> Option<Ciphertext<P>>
    where
        P: PointAdd,
    {
        if recipient_public.is_identity() {
            return None;
        }

        let ephemeral_public = base_point.scalar_mul(ephemeral_secret.as_limbs(), curve);
        let shared = recipient_public.scalar_mul(ephemeral_secret.as_limbs(), curve);
        let blinded_message = message.add(&shared, curve);

        Some(Ciphertext {
            ephemeral_public,
            blinded_message,
        })
    }

    /// Decrypt an ElGamal ciphertext back to the original point message.
    pub fn decrypt<P, const LIMBS: usize>(
        recipient_secret: &SecretScalar<LIMBS>,
        ciphertext: &Ciphertext<P>,
        curve: &P::Curve,
    ) -> P
    where
        P: PointAdd,
    {
        let shared = ciphertext
            .ephemeral_public
            .scalar_mul(recipient_secret.as_limbs(), curve);
        ciphertext.blinded_message.add(&shared.negate(curve), curve)
    }
}
