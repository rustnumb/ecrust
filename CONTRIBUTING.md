# Contributing to ecRust

Thank you for your interest in contributing to **ecRust**.

This project is a Rust workspace for experimenting with finite fields, elliptic curves, isogenies, and related cryptographic protocols. Contributions are welcome in code, tests, documentation, examples, benchmarks, and design discussions.

## Before you start

Please keep the following in mind:

- This project is **experimental**.
- Code in this repository should **not** be assumed to be production-ready or fully hardened against side-channel attacks unless explicitly documented otherwise.
- Correctness, clarity, and test coverage matter more than adding a large amount of code quickly.
- Small, focused pull requests are much easier to review than broad refactors.

## Ways to contribute

You can help by contributing:

- bug fixes
- new field or curve functionality
- tests and test vectors
- examples and demos
- documentation improvements
- benchmarks and profiling results
- API cleanup and ergonomics improvements

If you plan to make a large change, please open an issue or discussion first so the direction can be aligned early.

## Repository layout

The repository is organized as a Cargo workspace with multiple crates. The exact layout may evolve, but contributions will usually touch one or more of the following areas:

- `fp/` — finite field arithmetic
- `ec/` — elliptic-curve models and point arithmetic
- `isogeny/` — isogeny-related functionality
- `protocol/` — higher-level cryptographic protocols built on top of the algebraic layers

When possible, keep changes localized to the crate they belong to.

## Development setup

Clone the repository and build the whole workspace:

```bash
git clone https://github.com/rustnumb/ecrust
cd ecrust
cargo build --workspace
```

Run the full test suite:

```bash
cargo test --workspace
```

Generate documentation:

```bash
cargo doc --workspace --no-deps
```

Run formatting and lints before opening a pull request:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets --all-features -- -D warnings
```

If you are working on a single crate, you can limit commands accordingly, for example:

```bash
cargo test -p fp
cargo test -p ec
```

## Coding guidelines

### Rust style

- Format code with `cargo fmt`.
- Prefer simple, explicit code over clever code.
- Keep functions small when practical.
- Use descriptive names for algebraic objects and intermediates.
- Avoid unnecessary allocations in performance-sensitive paths.
- Keep public APIs consistent across crates where possible.

### Documentation

Public items should be documented.

Please add or update:

- `///` docs for public traits, structs, enums, and methods
- examples when an API is not obvious
- crate-level docs when adding a new module with public surface area

If a contribution changes behavior, also update any relevant:

- `README.md`
- examples
- demo files
- API docs

### Testing

Every non-trivial change should include tests.

Preferred kinds of tests include:

- deterministic unit tests
- cross-checks between equivalent formulas or models
- small-field exhaustive tests when feasible
- property-style tests for algebraic identities
- regression tests for previously failing cases

For arithmetic code, good tests often check:

- identity behavior
- inverses
- associativity on small samples
- doubling vs. self-addition
- scalar multiplication consistency
- on-curve validation after random generation

### Benchmarks and performance changes

Performance improvements are welcome, but correctness comes first.

If a pull request changes performance-sensitive code, it helps to include:

- a short explanation of the intended optimization
- benchmark numbers, if available
- any tradeoffs in readability, memory, or constant-time behavior

## Cryptographic and security-related contributions

This repository touches cryptographic code, so extra care is needed.

When contributing cryptographic functionality:

- do not assume an implementation is constant-time unless it has been reviewed with that goal in mind
- document any known limitations
- avoid introducing secret-dependent branches, memory accesses, or obvious timing leaks in code intended for hardened use
- separate experimental code from code meant to support stronger security claims
- explain assumptions clearly in docs and pull requests

If your change is only mathematically correct but not hardened, that is still useful — just say so explicitly.

## Pull request checklist

Before submitting a pull request, please check that:

- the code builds with `cargo build --workspace`
- tests pass with `cargo test --workspace`
- formatting is clean with `cargo fmt --all`
- clippy passes for the affected code
- public APIs are documented
- new behavior is covered by tests
- commit history is reasonably clean
- the pull request description explains what changed and why

A good pull request description usually includes:

- the problem being solved
- the approach taken
- any API changes
- any follow-up work left for later

## Commit messages

Please prefer clear commit messages.

Good examples:

- `fp: add tests for cubic extension inversion`
- `ec: implement Display for Jacobi quartic points`
- `isogeny: fix scalar edge case in demo`

## Issues

Bug reports and feature requests are welcome.

For bug reports, please include:

- what you expected
- what happened instead
- minimal reproduction steps
- relevant error messages
- the crate and module involved

For implementation questions, sketches, or design proposals, an issue or discussion is often the best starting point.

## License and contribution terms

By submitting a contribution, you agree that your contribution may be distributed under the same license terms as this repository.

If licensing expectations for the repository are updated or clarified later, this file should be kept in sync.

## Questions

If something is unclear, open an issue and ask. Early questions are much better than late rework.

Thank you for helping improve ecrust.
