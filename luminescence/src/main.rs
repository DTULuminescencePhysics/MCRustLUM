// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Command-line entry point for a luminescence Monte Carlo run.
//!
//! The program creates a directory under `run/`, copies `input.toml` into it,
//! and runs the simulation from that directory.

use std::error::Error;
use std::ffi::OsString;
use std::time::Instant;


/// Prepare the requested run directory and execute the default simulation workflow.
fn main() -> Result<(), Box<dyn Error>> {
    let folder_name = folder_name_from_arguments()?;
    io::filesystem::prepare_experiment_directory(folder_name.as_deref())?;
    common::random::set_seed(0);
    monte_carlo_run()
}

/// Load the copied input, run all repetitions, and consolidate their outputs.
///
/// The current executable requests ten repetitions, one parameter experiment,
/// and a minimum spatial ensemble of 25 traps.
fn monte_carlo_run() -> Result<(), Box<dyn Error>> {
    let total_start = Instant::now();
    let start = Instant::now();

    let inputs = io::inputs::read_inputs("input.toml")?;
    eprintln!("read inputs:       {:?}", start.elapsed());

    let start = Instant::now();
    let monte_carlo = mc::system_setup::MonteCarloSimulation::new(inputs, 10, 1,595)?;
    eprintln!("system setup:      {:?}", start.elapsed());
    
    let start = Instant::now();
    monte_carlo.run()?;
    eprintln!("Monte Carlo:       {:?}", start.elapsed());

    let start = Instant::now();
    mc::average::average_fill()?;
    eprintln!("average fill:      {:?}", start.elapsed());

    let start = Instant::now();
    mc::average::average_events()?;
    eprintln!("average events:    {:?}", start.elapsed());

    let start = Instant::now();

    let event_bin = Some(1.0);
    io::plotting::plot_default_results("average_fill.csv", "average_event.csv", ".", event_bin)?;
    eprintln!("plotting:          {:?}", start.elapsed());
    eprintln!("total:             {:?}", total_start.elapsed());

    Ok(())
}

/// Parse the optional single experiment-directory name from command-line arguments.
fn folder_name_from_arguments() -> Result<Option<OsString>, Box<dyn Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let folder_name = arguments.next();

    if arguments.next().is_some() {
        return Err("usage: luminescence [experiment-folder-name]".into());
    }

    Ok(folder_name)
}
