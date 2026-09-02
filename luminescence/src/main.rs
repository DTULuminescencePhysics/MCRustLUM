// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Command-line entry point for a luminescence Monte Carlo run.
//!
//! The program creates a directory under `run/`, copies `input.toml` into it,
//! and runs the simulation from that directory.

use std::error::Error;
use std::ffi::OsString;

fn main() -> Result<(), Box<dyn Error>> {
    monte_carlo_run()
}

fn monte_carlo_run() -> Result<(), Box<dyn Error>> {
    let folder_name = folder_name_from_arguments()?;
    io::filesystem::prepare_experiment_directory(folder_name.as_deref())?;

    let inputs = io::read_inputs("input.toml")?;
    let monte_carlo = mc::system_setup::MonteCarloSimulation::new(inputs, 10, 1)?;

    monte_carlo.run()?;

    Ok(())
}

fn folder_name_from_arguments() -> Result<Option<OsString>, Box<dyn Error>> {
    let mut arguments = std::env::args_os().skip(1);
    let folder_name = arguments.next();

    if arguments.next().is_some() {
        return Err("usage: luminescence [experiment-folder-name]".into());
    }

    Ok(folder_name)
}
