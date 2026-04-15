//! Generic isogeny abstraction



use fp::field_ops::FieldOps;
use ec::point_ops::PointOps;
use ec::curve_ops::Curve;
use subtle::{Choice, ConditionallySelectable};

pub trait IsogenyOps: Clone + ConditionallySelectable {
    type BaseField: FieldOps;
    type DomainCurve: Curve<BaseField = Self::BaseField>;
    type CodomainCurve: Curve<BaseField = Self::BaseField>;
    type DomainPoint: PointOps<Curve = Self::DomainCurve, BaseField = Self::BaseField>;
    type CodomainPoint: PointOps<Curve = Self::CodomainCurve, BaseField = Self::BaseField>;

    fn domain(&self) -> &Self::DomainCurve;
    fn codomain(&self) -> &Self::CodomainCurve;
    fn degree(&self) -> u64;
    fn evaluate(&self, p: &Self::DomainPoint) -> Self::CodomainPoint;
    fn is_separable(&self) -> Choice;
    //fn dual(&self) -> Option<Self>;
    //fn compose(&self, other: &Self) -> Option<Self>;
}

pub trait DualIsogenyOps: IsogenyOps {
    type Dual: IsogenyOps<
        BaseField = Self::BaseField,
        DomainCurve = Self::CodomainCurve,
        CodomainCurve = Self::DomainCurve,
        DomainPoint = Self::CodomainPoint,
        CodomainPoint = Self::DomainPoint,
    >;

    fn dual(&self) -> Self::Dual;
}