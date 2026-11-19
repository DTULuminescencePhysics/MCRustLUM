// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Infer thermal histories by comparing Monte Carlo predictions with measurements.
//!
//! [`thermo_chronology_experiment`] loads configuration, proposes time-temperature
//! histories, simulates final fill ratios, and stores accepted profiles.
//! Independent sampling and reversible-jump walks are supported; eligible
//! Bayesian walks can adapt proposal scales during burn-in. [`analyse`] builds
//! occupancy grids and representative paths for the final plots.
//!
//! For manual orchestration, construct [`system_setup::ChronologyExperiment`].

#![warn(missing_docs)]

/// Thermal-history records, constraints, and independent/walk proposal generation.
pub mod profiles;

/// Measurement residuals, statistical acceptance, and accepted-profile buffering.
pub mod acceptance;

/// Windowed acceptance monitoring and adaptive proposal-scale tuning.
pub mod burn_in;

/// Construction and execution of the inference workflow.
pub mod system_setup;
/// Analyses the results
pub mod analyse;

use std::{time::Instant, path::Path, error::Error};
use crate::analyse::{build_profile_grid,plot_profile_grid_heatmap};
use crate::profiles::TimeTempProfile;
use crate::system_setup::ChronologyExperiment;
use io::errors::PlotError;
use io::chrono_plotting::ProfileHeatmapOptions;
use io::inputs::{read_inputs, ResolvedExperiment};
use io::chrono_inputs::{ChronologyInputs,read_chrono_inputs};
use io::chrono_outputs::write_profile_grid_csv;

/// Run profile inference from `input.toml` and `chron_input.toml` in the current directory.
///
/// For each of `trial_num` proposals, simulate `experiments` configurations with
/// `repetitions` independent trajectories each. Write accepted profiles under
/// `tmp/`, then create grid and plot outputs. `use_real_profile` overlays the
/// first input profile. Configuration, simulation, and output errors propagate.
pub fn thermo_chronology_experiment(
    repetitions: usize, 
    experiments: usize, 
    trial_num: usize,
    use_real_profile: bool
)-> Result<(), Box<dyn Error>> 
{
    let start = Instant::now();

    let mc_inputs = read_inputs("input.toml", experiments)?;
    let chron_inputs = read_chrono_inputs("chron_input.toml")?;
    eprintln!("read inputs:       {:?}", start.elapsed());
    
    
    thermo_chronology_setup_run(&chron_inputs, &mc_inputs,repetitions,experiments,trial_num)?;
    let profile: Option<TimeTempProfile> = if use_real_profile {
        Some( TimeTempProfile {times: mc_inputs[0].time_temperature.times.value.clone(), temps: mc_inputs[0].time_temperature.temperatures.value.clone()  })
    } else{
       None
    };

    create_chronology_results(&chron_inputs,profile)?;

    
    Ok(())
}

/// Construct the inference state and run it with a deterministic derived seed.
///
/// Uses `repetitions * experiments * trial_num + 1` as the profile-generator
/// stream index and reports setup/run timing. Setup or simulation errors propagate.
fn thermo_chronology_setup_run(
    chron_inputs: &ChronologyInputs,
    mc_inputs: &Vec<ResolvedExperiment>,
    repetitions: usize, 
    experiments: usize, 
    trial_num: usize,
)-> Result<(), Box<dyn Error>> 
{
    
    let start = Instant::now();
    
    let mut chronology = ChronologyExperiment::setup(
        chron_inputs, mc_inputs, 
        repetitions, experiments)?;
    eprintln!("system setup:      {:?}", start.elapsed());
    
    let start = Instant::now();
    let seed = repetitions*experiments*trial_num +1;
    chronology.run(chron_inputs,trial_num, seed)?;
    eprintln!("Chronology:       {:?}", start.elapsed());

    Ok(())
}

/// Summarize stored accepted profiles into CSV and PNG files in the current directory.
///
/// Reads `tmp/chronology_results.bin.gz` and writes `profile_grid.csv`,
/// `accepted_profiles.png`, and `profile_heatmap.png`. An optional known
/// profile is overlaid in the stored units. Read, analysis, and plotting errors propagate.
pub fn create_chronology_results(
    profile_inputs: &ChronologyInputs, 
    real_profile: Option<TimeTempProfile>
) -> Result<(), Box<dyn Error>> 
{
    
    let output_csv = "profile_grid.csv";
    let image = "profile_heatmap.png";

    construct_heatmap_csv(profile_inputs, Path::new("tmp"), 
    "chronology_results.bin.gz", output_csv)?;
    
    let known_profile = real_profile.as_ref().map(|profile| {
        io::chrono_plotting::ProfilePath {
            times: &profile.times,
            temperatures: &profile.temps,
        }
    });

    io::chrono_plotting::plot_accepted_profiles::<TimeTempProfile>(
        Path::new("tmp").join("chronology_results.bin.gz"),
        "accepted_profiles.png",
        known_profile,
    )?;

    plot_heatmap(real_profile, output_csv, image)?;
    
    Ok(())

}

/// Read accepted-profile batches, build an occupancy grid, and write it to CSV.
///
/// The source is `output_directory/file_name`; `output_csv` is used as supplied,
/// relative to the current directory when not absolute. Grid settings come
/// from `profile_inputs`. Invalid data and file errors propagate.
pub fn construct_heatmap_csv(
    profile_inputs: &ChronologyInputs,
    output_directory: impl AsRef<Path>, 
    file_name: &str, output_csv: &str
) -> Result<(), Box<dyn Error>> 
{

    let output_directory: &Path = output_directory.as_ref();
    let output_file = output_directory.join(file_name);

    let grid = build_profile_grid(output_file,
        &profile_inputs.grid_options)?;
    
    
    let results_file = Path::new(output_csv);
    
    write_profile_grid_csv(results_file, &grid)?;
    Ok(())

}


/// Render a profile-grid CSV using default heatmap settings.
///
/// An optional known profile must use the CSV units. Writes to `file_name`
/// and returns errors for invalid data or drawing/output failures.
pub fn plot_heatmap(real_profile: Option<TimeTempProfile>, csv_path: impl AsRef<Path>, file_name: &str) -> Result< (), PlotError> {

    let output_path = Path::new(file_name);
    let options = ProfileHeatmapOptions::default();
    
    match real_profile {
        Some(real_profile) => {
            plot_profile_grid_heatmap(
                                    csv_path, 
                                    output_path, 
                                    &options, 
                                    Some(&real_profile))
        },
        None => {
            plot_profile_grid_heatmap(
                csv_path, 
                output_path, 
                &options, 
                None)
        }
    }
}




