// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Command-line entry point for a luminescence Monte Carlo run.
//!
//! The program reads `input.toml` from the working directory when present and
//! otherwise constructs the simulation from the library defaults.

use std::{error::Error, path::Path};

fn main() -> Result<(), Box<dyn Error>> {
    monte_carlo_run()
}

fn monte_carlo_run() -> Result<(), Box<dyn Error>> {
    let inputs = if Path::new("input.toml").exists() {
        io::read_inputs("input.toml")?
    } else {
        io::default_inputs()
    };
    let monte_carlo = mc::system_setup::MonteCarloSimulation::new(inputs, 10, 1)?;

    monte_carlo.run()?;

    Ok(())
}
