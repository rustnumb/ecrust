use crypto_bigint::{Uint, const_prime_monty_params};

use ec::curve_weierstrass::WeierstrassCurve;
use ec::point_weierstrass::AffinePoint;
use fp::fp_element::FpElement;
use protocol::ecdh::Ecdh;
use protocol::elgamal::EcElGamal;
use protocol::scalar::SecretScalar;

// SEC 2 secp521r1 domain parameters.
const_prime_monty_params!(
    Secp521r1Mod,
    Uint<9>,
    "00000000000001FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
    1
);
type Fp521 = FpElement<Secp521r1Mod, 9>;

fn curve() -> WeierstrassCurve<Fp521> {
    let a = Fp521::from_words([
        0xfffffffffffffffc,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x1ff,
    ]);

    let b = Fp521::from_words([
        0xef451fd46b503f00,
        0x3573df883d2c34f1,
        0x1652c0bd3bb1bf07,
        0x56193951ec7e937b,
        0xb8b489918ef109e1,
        0xa2da725b99b315f3,
        0x929a21a0b68540ee,
        0x953eb9618e1c9a1f,
        0x0000000000000051,
    ]);

    WeierstrassCurve::new_short(a, b)
}

fn base_point() -> AffinePoint<Fp521> {
    let x = Fp521::from_words([
        0xf97e7e31c2e5bd66,
        0x3348b3c1856a429b,
        0xfe1dc127a2ffa8de,
        0xa14b5e77efe75928,
        0xf828af606b4d3dba,
        0x9c648139053fb521,
        0x9e3ecb662395b442,
        0x858e06b70404e9cd,
        0x00000000000000c6,
    ]);

    let y = Fp521::from_words([
        0x88be94769fd16650,
        0x353c7086a272c240,
        0xc550b9013fad0761,
        0x97ee72995ef42640,
        0x17afbd17273e662c,
        0x98f54449579b4468,
        0x5c8a5fb42c7d1bd9,
        0x39296a789a3bc004,
        0x0000000000000118,
    ]);

    AffinePoint::new(x, y)
}

#[test]
fn ecdh_shared_secret_matches() {
    let curve = curve();
    let base = base_point();

    let alice_sk = SecretScalar::<9>::new([5, 0, 0, 0, 0, 0, 0, 0, 0]);
    let bob_sk = SecretScalar::<9>::new([7, 0, 0, 0, 0, 0, 0, 0, 0]);

    let alice_pk = Ecdh::derive_public_key(&base, &alice_sk, &curve);
    let bob_pk = Ecdh::derive_public_key(&base, &bob_sk, &curve);

    let alice_ss = Ecdh::shared_secret(&alice_sk, &bob_pk, &curve)
        .expect("peer public key must be non-identity");
    let bob_ss = Ecdh::shared_secret(&bob_sk, &alice_pk, &curve)
        .expect("peer public key must be non-identity");

    assert_eq!(alice_ss, bob_ss);
}

#[test]
fn elgamal_roundtrip_on_curve_point() {
    let curve = curve();
    let base = base_point();

    let recipient_sk = SecretScalar::<9>::new([9, 0, 1, 0, 1, 0, 0, 1, 0]);
    let recipient_pk = Ecdh::derive_public_key(&base, &recipient_sk, &curve);

    let message = base.scalar_mul(&[4, 0, 0, 0, 0, 1, 0, 0, 0], &curve);
    let nonce = SecretScalar::<9>::new([3, 0, 0, 1, 0, 0, 0, 0, 0]);

    let ciphertext = EcElGamal::encrypt(&base, &recipient_pk, &message, &nonce, &curve)
        .expect("recipient public key must be non-identity");
    let decrypted = EcElGamal::decrypt(&recipient_sk, &ciphertext, &curve);

    assert_eq!(decrypted, message);
}
