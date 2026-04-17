//! Umbrella crate for the `ecrust` ecosystem.
//!
//! This crate re-exports the modular workspace crates behind one package.

#[cfg(feature = "fp")]
pub use fp;

#[cfg(feature = "ec")]
pub use ec;

#[cfg(feature = "isogeny")]
pub use isogeny;

#[cfg(feature = "protocol")]
pub use protocol;
