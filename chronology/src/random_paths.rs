// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use common::numeric::Float;
use crate::profiles::{TimeTempProfile,ProfileType};
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use mc::system_setup::MonteCarloSimulation;
use io::chrono_inputs::RandomWalkInputs;
use io::inputs::SimulationInputs;
use io::chrono_outputs;
use rand::{Rng, RngExt};
use rand::distr::Open01;
use std::path::Path;


pub struct RandomWalkSimulation{

    fill_ratio: Vec<Float>,
    fill_ratio_error: Vec<Float>, 
    profile_type: ProfileType,
    monte_carlo:MonteCarloSimulation,
    time_unit:TimeUnit,
    temp_unit:TemperatureUnit,
}

impl  RandomWalkSimulation {

    pub fn new(fill_ratio: Vec<Float>, 
        fill_ratio_error: Vec<Float>, 
        profile_inputs: RandomWalkInputs,
        inputs: &SimulationInputs,
        repetitions: usize,
        experiments: usize,
        minimum_traps:usize,
    ) -> Result<Self, String> {

        let profile_type = ProfileType::new_random_walk(&profile_inputs)?;
        let monte_carlo = MonteCarloSimulation::new(inputs.clone(), repetitions, experiments, minimum_traps)?;
        
        Ok(
            Self {
                fill_ratio,
                fill_ratio_error,
                profile_type,
                monte_carlo,
                time_unit: profile_inputs.time_unit,
                temp_unit: profile_inputs.temp_unit,
            }
        )
}

    pub fn run_random_walks(&self, trial_num: usize, output_directory: impl AsRef<Path>, rng: &mut impl Rng) -> Result<(),String> {
        let output_directory = output_directory.as_ref();
        let output_path = output_directory.join(format!(
                "chronology_results.bin.gz"
            ));
        
        chrono_outputs::create_temporary_chronology_experiment_file(&output_path)
                                                                    .map_err(|error| error.to_string())?;
           
        chrono_outputs::write_temporary_chronology_units(&output_path, self.time_unit, self.temp_unit)
                                                        .map_err(|error| error.to_string())?;
        
        let mut accepted_profile: Vec<TimeTempProfile> = Vec::with_capacity(1000);
        
        for num in 0.. trial_num {
            let tt_prof = self.profile_type.create_new_profile(rng)?;
            let result = self.monte_carlo.run_to_final_ratio_only(tt_prof.generate_profile(self.time_unit,self.temp_unit)?, num)?;
            if self.check_profile(result, rng){
                accepted_profile.push(tt_prof);
            }
            if accepted_profile.len() == 1000 {
                chrono_outputs::append_chronology_profile_to_experiment_file( &output_path, &accepted_profile)
                                                                            .map_err(|error| error.to_string())?;
                accepted_profile.clear(); 
            }
        }
        if accepted_profile.len() != 0 {
            chrono_outputs::append_chronology_profile_to_experiment_file( &output_path, &accepted_profile)
                                                                            .map_err(|error| error.to_string())?;
                accepted_profile.clear();
        }

        Ok(())
    }

    pub fn check_profile(&self, result: Vec<Vec<Float>>, rng: &mut impl Rng) -> bool{
        let mut total: Float = 0.0;

        for ((res, target), error) in result.iter()
                                                                   .zip(&self.fill_ratio)
                                                                   .zip(&self.fill_ratio_error)
        {
            let diff: Float = ((res.iter().sum::<Float>())/(res.len()as Float-target))/error;
            total += diff.powf(2.0);
        }

        if (-0.5*total).exp() > rng.sample(Open01){
            return true; 
        }else{
            return false; 
        }
    }

    
}






