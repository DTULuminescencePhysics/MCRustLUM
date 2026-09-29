// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Loading and representing simulation configuration.
//!
//! Use [`inputs::read_inputs`] when a TOML file is supplied and
//! [`inputs::default_inputs`]
//! when the built-in configuration is sufficient. Both paths return the same
//! [`inputs::SimulationInputs`] type, so downstream simulation code does not need to
//! know where the values came from.


pub mod errors;

/// Typed groups corresponding to the sections of an input TOML file.
pub mod inputs;

/// Creation of per-run experiment and temporary-output directories.
pub mod filesystem;

/// Writing consolidated, user-facing simulation output.
pub mod outputs;

/// Loading consolidated CSV results and producing plots from them.
pub mod plotting;


pub mod chrono_inputs;

pub mod chrono_outputs;


