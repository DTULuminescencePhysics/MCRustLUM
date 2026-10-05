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

    let total_start = Instant::now();

    let start = Instant::now();
    mc::monte_carlo_run(10, 1, 400, None)?;
    eprintln!("Monte Carlo Complete:             {:?}", start.elapsed());
    
    // let start = Instant::now();
    // chronology::thermo_chronology_run(10,1,595,100,false)?;
    // eprintln!("Thermochronology Complete:             {:?}", start.elapsed());
    
    
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
