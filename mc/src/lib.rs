// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only 

//! Kinetic Monte Carlo model setup and event-time generation.
//!
//! [`system_setup::MonteCarloSimulation`] turns the grouped input values from
//! the `io` crate into the crystal, time/temperature profile, and transition
//! selections used during a run.

/// Construction and resetting of Monte Carlo simulation state.
pub mod system_setup;
/// Contains experiment specific code
pub mod experiment;


