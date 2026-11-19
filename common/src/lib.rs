// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Shared simulation types.
//!
//! This crate contains, unit conversions the typed inputs used to evaluate the rate
//! equations and output parameters. File parsing and Monte Carlo orchestration live in the sibling
//! `io` and `mc` crates while specific rate equation is in `physics`.

#![warn(missing_docs)]

/// Numeric precision aliases and shape-aware arithmetic traits.
pub mod numeric;
/// Physical constants and metric, time, temperature conversions.
pub mod constants;
/// Holds unit structs and helper routines needed to make conversions.
pub mod units;
/// Holds structs describing specific quantities and the methods for how they can be generated.
pub mod quantities;
/// Compact identifiers and occupancy partitions for traps, holes, and band tails.
pub mod place_ids;
/// Deterministic base-seed handling and independent repetition generators.
pub mod random;

