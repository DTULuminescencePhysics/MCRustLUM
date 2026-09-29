// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

pub mod profiles;

pub mod random_paths;



use random_paths::RandomWalkSimulation;
use io::chrono_inputs::ChronologyInputs;
use io::inputs::SimulationInputs;
use rand::Rng;
use std::path::Path;
pub enum ChronologyExperiment{

    RandomWalk {experiment: RandomWalkSimulation },

    RJMCMC,

}

impl ChronologyExperiment {

    pub fn setup( profile_inputs:ChronologyInputs, 
                  inputs: SimulationInputs, 
                  repetitions: usize, 
                  experiments: usize, 
                  minimum_traps:usize) 
        -> Result<Self,String>
    {

        if profile_inputs.setup.random_paths {
            let experiment = RandomWalkSimulation::new(profile_inputs.setup.fill_ratio, 
                                                                            profile_inputs.setup.fill_ratio_error, 
                                                                            profile_inputs.random_walks, 
                                                                            inputs, 
                                                                            repetitions,
                                                                            experiments, 
                                                                            minimum_traps)?;
            return Ok(Self::RandomWalk { experiment });
        }else{
            return Err(format!("RJMCM not set up yet"))
        }

    }

    pub fn run(&self, trial_num: usize, rng: &mut impl Rng) -> Result<(), String> {
        
        match self {
            ChronologyExperiment::RandomWalk { experiment } => {
                experiment.run_random_walks(trial_num, Path::new("tmp"), rng)?;
                Ok(())
            },
            ChronologyExperiment::RJMCMC => {
                return Err(format!("RJMCM not set up yet"))

            }
        }

    }


}