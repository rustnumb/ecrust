//! Generic isogeny abstraction.
//!
//! # Overview
//!
//! An isogeny is a morphism of elliptic curves
//!
//! ```text
//! φ : E -> E'
//! ```
//!
//! that is also a group homomorphism.
//!
//! This trait is meant to capture the **minimal common interface**
//! that any concrete isogeny implementation should expose:
//!
//! - its base field,
//! - its domain and codomain curves,
//! - the point representations used on those curves,
//! - its degree,
//! - evaluation on points,
//! - and whether it is separable.
//!
//! # Why separate `DomainPoint` and `CodomainPoint`?
//!
//! In principle, the source and target curves may use different native
//! point representations.
//!
//! For example:
//! - the domain could use affine Weierstrass points,
//! - the codomain could use projective points,
//! - or an x-only Montgomery/Kummer representation.
//!
//! So we keep both associated types explicit instead of assuming that the
//! same point type is used on both sides.

use ec::curve_ops::Curve;
use ec::point_ops::PointOps;
use fp::field_ops::FieldOps;
use subtle::Choice;

/// Generic interface for an isogeny.
///
/// An implementation of this trait represents a map
///
/// ```text
/// φ : E -> E'
/// ```
///
/// together with the information needed to evaluate it on points.
///
/// The trait is intentionally small: it only describes the common
/// functionality of isogenies, without forcing a particular internal
/// representation (kernel-based, Vélu formulas, x-only formulas, chains, etc.).
pub trait IsogenyOps: Clone {
    /// Base field over which both the domain and codomain curves are defined.
    type BaseField: FieldOps;

    /// Domain curve `E` of the isogeny `φ : E -> E'`.
    type DomainCurve: Curve<BaseField = Self::BaseField>;

    /// Codomain curve `E'` of the isogeny `φ : E -> E'`.
    type CodomainCurve: Curve<BaseField = Self::BaseField>;

    /// Native point representation used when evaluating points on the domain.
    type DomainPoint: PointOps<Curve = Self::DomainCurve, BaseField = Self::BaseField>;

    /// Native point representation produced when evaluating the isogeny.
    type CodomainPoint: PointOps<Curve = Self::CodomainCurve, BaseField = Self::BaseField>;

    /// Return the domain curve of the isogeny.
    fn domain(&self) -> &Self::DomainCurve;

    /// Return the codomain curve of the isogeny.
    fn codomain(&self) -> &Self::CodomainCurve;

    /// Return the degree of the isogeny.
    ///
    /// For a composition `ψ ∘ φ`, the degree should be
    ///
    /// ```text
    /// deg(ψ ∘ φ) = deg(ψ) deg(φ).
    /// ```
    fn degree(&self) -> u64;

    /// Evaluate the isogeny at a point of the domain curve.
    ///
    /// If `φ : E -> E'` and `P ∈ E`, this returns `φ(P) ∈ E'`.
    fn evaluate(&self, p: &Self::DomainPoint) -> Self::CodomainPoint;

    /// Return whether the isogeny is separable.
    ///
    /// In many protocol settings one mainly works with separable isogenies,
    /// but we keep this as an explicit property of the map.
    fn is_separable(&self) -> Choice;
}

/// Extension trait for isogenies whose dual isogeny is available.
///
/// If
///
/// ```text
/// φ : E -> E'
/// ```
///
/// is an isogeny, then its dual
///
/// ```text
/// φ^∨ : E' -> E
/// ```
///
/// goes in the opposite direction, has the same degree, and satisfies
/// the standard duality relations.
pub trait DualIsogenyOps: IsogenyOps {
    /// Type of the dual isogeny.
    ///
    /// Note how domain and codomain are swapped, and likewise the point
    /// representations on the source and target.
    type Dual: IsogenyOps<
            BaseField = Self::BaseField,
            DomainCurve = Self::CodomainCurve,
            CodomainCurve = Self::DomainCurve,
            DomainPoint = Self::CodomainPoint,
            CodomainPoint = Self::DomainPoint,
        >;

    /// Return the dual isogeny.
    fn dual(&self) -> Self::Dual;
}

/// A wrapper representing the composition of two isogenies.
///
/// If
///
/// ```text
/// first  : E  -> E'
/// second : E' -> E'',
/// ```
///
/// then `CompositeIsogeny { first, second }` represents
///
/// ```text
/// second ∘ first : E -> E''.
/// ```
///
/// # Why a wrapper?
///
/// This is a very convenient generic representation of composition:
/// instead of trying to "simplify" the composite into a single special
/// formula, we just store the two maps and evaluate them one after the other.
///
/// This is often the cleanest first design, especially when different
/// concrete isogeny families may be composed together.
#[derive(Clone)]
pub struct CompositeIsogeny<I1, I2> {
    /// The first isogeny in the composition.
    pub first: I1,

    /// The second isogeny in the composition.
    pub second: I2,
}

impl<I1, I2> IsogenyOps for CompositeIsogeny<I1, I2>
where
    // The first map is an arbitrary isogeny
    I1: IsogenyOps,
    // The second map must start where the first one ends.
    //
    // Concretely, if
    //   first  : E  -> E'
    // then we require
    //   second : E' -> E''
    //
    // This means:
    // - same base field,
    // - domain curve of `second` = codomain curve of `first`,
    // - domain point type of `second` = codomain point type of `first`.
    I2: IsogenyOps<
            BaseField = I1::BaseField,
            DomainCurve = I1::CodomainCurve,
            DomainPoint = I1::CodomainPoint,
        >,
{
    /// The composite map is still defined over the same base field.
    type BaseField = I1::BaseField;

    /// The domain of `second ∘ first` is the domain of `first`.
    type DomainCurve = I1::DomainCurve;

    /// The codomain of `second ∘ first` is the codomain of `second`.
    type CodomainCurve = I2::CodomainCurve;

    /// Points fed into the composite are points of the domain of `first`.
    type DomainPoint = I1::DomainPoint;

    /// Points output by the composite are points of the codomain of `second`.
    type CodomainPoint = I2::CodomainPoint;

    fn domain(&self) -> &Self::DomainCurve {
        self.first.domain()
    }

    fn codomain(&self) -> &Self::CodomainCurve {
        self.second.codomain()
    }

    fn degree(&self) -> u64 {
        // By multiplicativity of degree:
        // deg(second ∘ first) = deg(second) * deg(first)
        self.first.degree() * self.second.degree()
    }

    fn evaluate(&self, p: &Self::DomainPoint) -> Self::CodomainPoint {
        // First apply `first : E -> E'`
        let q = self.first.evaluate(p);

        // Then apply `second : E' -> E''`
        self.second.evaluate(&q)
    }

    fn is_separable(&self) -> Choice {
        // A composition of separable isogenies is separable.
        // So the composite is separable exactly when both factors are.
        self.first.is_separable() & self.second.is_separable()
    }
}
