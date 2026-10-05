// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

pub mod profiles;

pub mod random_paths;

pub mod system_setup;

pub mod analyse;

use std::time::Instant;
use crate::analyse::{build_profile_grid,plot_profile_grid_heatmap};
use crate::profiles::TimeTempProfile;
use crate::system_setup::ChronologyExperiment;
use std::error::Error;
use io::errors::PlotError;
use io::chrono_plotting::ProfileHeatmapOptions;
use io::inputs::{SimulationInputs, read_inputs};
use io::chrono_inputs::{ChronologyInputs,read_chrono_inputs};
use io::chrono_outputs::write_profile_grid_csv;
use std::path::Path;

fn thermo_chronology_setup(
    chron_inputs: &ChronologyInputs,
    mc_inputs: &SimulationInputs,
    repetitions: usize, 
    experiments: usize, 
    minimum_traps: usize,
    trial_num: usize,
)-> Result<(), Box<dyn Error>> 
{
    
    let start = Instant::now();
    
    let chronology = ChronologyExperiment::setup(
        chron_inputs, mc_inputs, 
        repetitions, experiments, minimum_traps)?;
    eprintln!("system setup:      {:?}", start.elapsed());
    
    let start = Instant::now();
    let seed = repetitions*experiments*trial_num +1;
    chronology.run(trial_num, seed)?;
    eprintln!("Chronology:       {:?}", start.elapsed());

    Ok(())
}

pub fn thermo_chronology_run(
    repetitions: usize, 
    experiments: usize, 
    minimum_traps: usize,
    trial_num: usize,
    use_real_profile: bool
)-> Result<(), Box<dyn Error>> 
{
    let start = Instant::now();

    let mc_inputs = read_inputs("input.toml")?;
    let chron_inputs = read_chrono_inputs("chron_input.toml")?;
    eprintln!("read inputs:       {:?}", start.elapsed());
    
    
    thermo_chronology_setup(&chron_inputs, &mc_inputs,repetitions,experiments,minimum_traps,trial_num)?;
    chronology_result(&chron_inputs,use_real_profile)?;

    
    Ok(())
}

pub fn chronology_result(chron_inputs: &ChronologyInputs, use_real_profile: bool) -> Result<(), Box<dyn Error>>{
    
    let profile = if use_real_profile {
        None
    } else{
        None
    };
    create_chronology_results(chron_inputs,profile)?;

    Ok(())

}

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

pub fn plot_heatmap(real_profile: Option<&TimeTempProfile>, csv_path: impl AsRef<Path>, file_name: &str) -> Result< (), PlotError> {

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

pub fn create_chronology_results(
    profile_inputs: &ChronologyInputs, 
    real_profile: Option<&TimeTempProfile>
) -> Result<(), Box<dyn Error>> 
{
    
    let output_csv = "profile_grid.csv";
    let image = "profile_heatmap.png";

    construct_heatmap_csv(profile_inputs, Path::new("tmp"), 
    "chronology_results.bin.gz", output_csv)?;
    
    plot_heatmap(real_profile, output_csv, image)?;
    
    Ok(())


}