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
use common::numeric::TimeFloat;

use io::inputs::{read_inputs, SimulationInputs};
use io::plotting;

/// Main call function to run the required Monte Carlo experiments

pub fn monte_carlo_run(repetitions: usize, experiments: usize, minimum_traps: usize, event_bin_width: Option<TimeFloat>) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();

    let inputs = read_inputs("input.toml")?;
    eprintln!("read inputs:       {:?}", start.elapsed());

    monte_carlo_experiment(&inputs, repetitions, experiments, minimum_traps)?;

    monte_carlo_result(&inputs, experiments, event_bin_width)?;
    Ok(())

}

/// Main call function to run the required Monte Carlo experiments
fn monte_carlo_experiment(inputs: &SimulationInputs, repetitions: usize, experiments: usize, minimum_traps: usize) -> Result<(), Box<dyn Error>> {


    let start = Instant::now();
    let monte_carlo = MonteCarloSimulation::new(inputs.clone(), repetitions, experiments, minimum_traps)?;
    eprintln!("system setup:      {:?}", start.elapsed());

    let start = Instant::now();
    monte_carlo.run()?;
    eprintln!("Monte Carlo:       {:?}", start.elapsed());

    Ok(())
}

/// Function that creates the default set of results for the Monte Carlo experiments
fn monte_carlo_result(inputs: &SimulationInputs, experiments: usize, event_bin_width: Option<TimeFloat>) -> Result<(), Box<dyn Error>>{
    
    let start = Instant::now();
    let output_directory = ".";
    
    plotting::plot_time_temperature(
        output_directory,
        &inputs.time_temperature)?;


    crate::average::average_fill()?;
    eprintln!("average fill:      {:?}", start.elapsed());
    
    let start = Instant::now();
    crate::average::average_events()?;
    eprintln!("average events:    {:?}", start.elapsed());

    let start = Instant::now();

    for exp in 0..experiments{
        let fill_csv = format!("average_fill_{}.csv", exp);
        let event_csv = format!("average_event_{}.csv", exp); 
        
        plotting::plot_default_continuous_results(
            fill_csv,
            output_directory,
            inputs.time_temperature.time_unit,
            inputs.time_temperature.temp_unit,)?;
        
        plotting::plot_default_event_results(
            event_csv,
            output_directory,
            event_bin_width,
            inputs.time_temperature.time_unit,
        )?;
    }   
    
    eprintln!("plotting:          {:?}", start.elapsed());

    Ok(())
}
