// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Configuration loading, compressed result storage, CSV export, and plotting.
//!
//! Start with [`inputs::read_inputs`] to resolve named experiment definitions
//! and normalize their units into [`inputs::ResolvedExperiment`] values.
//! [`chrono_inputs::read_chrono_inputs`] loads inference settings separately.
//! Simulation trajectories are streamed through [`outputs`]; accepted thermal
//! histories use [`chrono_outputs`]. The plotting modules consume consolidated
//! results and report path-aware errors.

#![warn(missing_docs)]

/// Path-aware input, output, filesystem, and plotting error types.
pub mod errors;

/// Typed groups corresponding to the sections of an input TOML file.
pub mod inputs;

/// Creation of per-run experiment and temporary-output directories.
pub mod filesystem;

/// Writing consolidated, user-facing simulation output.
pub mod outputs;

/// Loading consolidated CSV results and producing plots from them.
pub mod plotting;

/// TOML settings for thermal-history inference and adaptive burn-in.
pub mod chrono_inputs;

/// Accepted-profile batch storage and occupancy-grid CSV export.
pub mod chrono_outputs;

/// Plot weighted time-temperature paths produced by chronology analysis.
pub mod chrono_plotting;
