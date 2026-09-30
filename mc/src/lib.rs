// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Kinetic Monte Carlo model setup and event-time generation.
//!
//! [`system_setup::MonteCarloSimulation`] turns the grouped input values from
//! the `io` crate into the crystal, time/temperature profile, and transition
//! selections used during a run.

#![warn(missing_docs)]

/// Per-repetition Monte Carlo state and experiment dispatch.
pub mod experiment;
/// Construction and resetting of Monte Carlo simulation state.
pub mod system_setup;

/// Event-rate construction, lifetime sampling, and kinetic Monte Carlo stepping.
pub mod calculate_times;

/// Streaming consolidation of repeated fill trajectories and event counts.
pub mod average;

use std::time::Instant;
use std::error::Error;

use crate::system_setup::MonteCarloSimulation;
use io::inputs::read_inputs;
use io::plotting;

/// Main call function to run the required Monte Carlo experiments
pub fn monte_carlo_run(repetitions: usize, experiments: usize, minimum_traps: usize) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();

    let inputs = read_inputs("input.toml")?;
    eprintln!("read inputs:       {:?}", start.elapsed());

    let start = Instant::now();
    let monte_carlo = MonteCarloSimulation::new(inputs, repetitions, experiments, minimum_traps)?;
    eprintln!("system setup:      {:?}", start.elapsed());

    let start = Instant::now();
    monte_carlo.run()?;
    eprintln!("Monte Carlo:       {:?}", start.elapsed());

    Ok(())
}

/// Function that creates the default set of results for the Monte Carlo experiments
pub fn monte_carlo_result() -> Result<(), Box<dyn Error>>{
    let start = Instant::now();
    crate::average::average_fill()?;
    eprintln!("average fill:      {:?}", start.elapsed());

    let start = Instant::now();
    crate::average::average_events()?;
    eprintln!("average events:    {:?}", start.elapsed());

    let start = Instant::now();

    let event_bin = Some(1.0);
    plotting::plot_default_results(
        "average_fill_0.csv",
        "average_event_0.csv",
        ".",
        event_bin,
    )?;
    eprintln!("plotting:          {:?}", start.elapsed());

    Ok(())
}
