use crypto_bigint::{Uint, const_prime_monty_params};

use ec::curve_edwards::EdwardsCurve;
use ec::curve_jacobi_intersection::JacobiIntersectionCurve;
use ec::curve_jacobi_quartic::JacobiQuarticCurve;
use ec::curve_montgomery::MontgomeryCurve;
use ec::curve_ops::Curve;
use ec::curve_weierstrass::WeierstrassCurve;

use fp::fp_element::FpElement;

use rand::rngs::ThreadRng;

// Small demo field: F_19
const_prime_monty_params!(Fp19Mod, Uint<1>, "0000000000000013", 2);
type F19 = FpElement<Fp19Mod, 1>;

fn fp(x: u64) -> F19 {
    F19::from_u64(x)
}

fn show_curve<C>(name: &str, curve: &C, rng: &mut ThreadRng)
where
    C: Curve + core::fmt::Display,
    C::Point: core::fmt::Display,
{
    println!("============================================================");
    println!("{name}");
    println!("------------------------------------------------------------");
    println!("curve compact : {}", curve);
    println!("curve pretty  :\n{:#}", curve);

    let p = curve.random_point(rng);

    println!("point compact : {}", p);
    println!("point pretty  :\n{:#}", p);
    println!("on curve?     : {}", curve.is_on_curve(&p));
    println!();
}

fn main() {
    let mut rng = rand::rng();
    // 1. Short Weierstrass over F_19: y^2 = x^3 + 2x + 3
    let w = WeierstrassCurve::new_short(fp(2), fp(3));
    show_curve("Weierstrass", &w, &mut rng);

    // 2. Montgomery over F_19: B y^2 = x(x^2 + A x + 1)
    // Smooth if B != 0 and A != ±2 in odd characteristic.
    let m = MontgomeryCurve::new(fp(3), fp(1));
    show_curve("Montgomery", &m, &mut rng);

    // 3. Edwards over F_19: x^2 + y^2 = 1 + d x^2 y^2
    // Pick d = 2 (nonzero, not 1; also a nonsquare in F_19).
    let e = EdwardsCurve::new(fp(2));
    show_curve("Edwards", &e, &mut rng);

    // 4. Jacobi quartic over F_19: y^2 = d x^4 + 2 a x^2 + 1
    // Need d != 0 and a^2 != d.
    let jq = JacobiQuarticCurve::new(fp(3), fp(5));
    show_curve("Jacobi quartic", &jq, &mut rng);

    // 5. Jacobi intersection over F_19:
    // s^2 + c^2 = 1,  a s^2 + d^2 = 1
    // Need a != 0, 1.
    let ji = JacobiIntersectionCurve::new(fp(2));
    show_curve("Jacobi intersection", &ji, &mut rng);
}
