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
use crate::average::{average_fill,rolling_average_fill,average_events};
use common::numeric::TimeFloat;

use io::inputs::{read_inputs, ResolvedExperiment};
use io::plotting::{plot_time_temperature, plot_default_continuous_results,plot_rolling_average_fill_results, plot_default_event_results};

/// Main call function to run the required Monte Carlo experiments

pub fn monte_carlo_run(repetitions: usize, experiments: usize, minimum_traps: usize, event_bin_width: Option<TimeFloat>) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();

    let inputs = read_inputs("input.toml", experiments)?;
    eprintln!("read inputs:       {:?}", start.elapsed());

    monte_carlo_experiment(&inputs, repetitions, experiments, minimum_traps)?;

    monte_carlo_result(&inputs, experiments, event_bin_width)?;
    Ok(())

}

/// Main call function to run the required Monte Carlo experiments
fn monte_carlo_experiment(inputs: &Vec<ResolvedExperiment>, repetitions: usize, experiments: usize, minimum_traps: usize) -> Result<(), Box<dyn Error>> {


    let start = Instant::now();
    let mut monte_carlo = MonteCarloSimulation::new(&inputs, repetitions, experiments, minimum_traps)?;
    eprintln!("system setup:      {:?}", start.elapsed());

    let start = Instant::now();
    monte_carlo.run()?;
    eprintln!("Monte Carlo:       {:?}", start.elapsed());

    Ok(())
}

/// Function that creates the default set of results for the Monte Carlo experiments
fn monte_carlo_result(inputs: &Vec<ResolvedExperiment>, experiments: usize, event_bin_width: Option<TimeFloat>) -> Result<(), Box<dyn Error>>{
    
    let start = Instant::now();
    let output_directory = ".";
    
    for input in inputs {
        plot_time_temperature(
        output_directory,
        &input.time_temperature
    )?;}
    

    average_fill()?;
    eprintln!("Average fill calculation:      {:?}", start.elapsed());
    
    let start = Instant::now();
    let rolling_divisor = 10;
    rolling_average_fill(rolling_divisor)?;
    eprintln!("Rolling average calculation:      {:?}", start.elapsed());
    
    let start = Instant::now();

    average_events()?;
    eprintln!("Average events calculation:    {:?}", start.elapsed());

    let start = Instant::now();

    for exp in 0..experiments{
        let fill_csv = format!("average_fill_{}.csv", exp);
        let av_fill_csv = format!("rolling_average_fill_{}.csv", exp);
        let event_csv = format!("average_event_{}.csv", exp); 
        
        let time_unit = inputs[exp].time_temperature.times.unit.get_time_unit()?;
        let temp_unit = inputs[exp].time_temperature.temperatures.unit.get_temperature_unit()?;
        plot_default_continuous_results(
            exp,
            fill_csv,
            output_directory,
            time_unit,           
            temp_unit,
        )?;

        plot_rolling_average_fill_results(
            av_fill_csv,
            output_directory,
            time_unit,

        )?;
        
        plot_default_event_results(
            event_csv,
            output_directory,
            event_bin_width,
            time_unit,
        )?;
    }   
    
    eprintln!("Plotting:          {:?}", start.elapsed());

    Ok(())
}
