// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Command-line entry point for a luminescence Monte Carlo run.
//!
//! The program creates a directory under `run/`, copies `input.toml` into it,
//! and runs the simulation from that directory.

use std::error::Error;
use std::ffi::OsString;

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
    let inputs = io::inputs::read_inputs("input.toml")?;
    let monte_carlo = mc::system_setup::MonteCarloSimulation::new(inputs, 40, 1,25)?;

    monte_carlo.run()?;
    mc::average::average_fill()?;
    mc::average::average_events()?;
    let event_bin = Some(1.0);
    io::plotting::plot_default_results("average_fill.csv", "average_event.csv", ".", event_bin)?;

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
