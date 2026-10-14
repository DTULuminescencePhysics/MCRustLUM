// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use crate::random_paths::RandomWalkSimulation;
use crate::markov_chain_mc::RJMCMCSimulation;
use io::chrono_inputs::ChronologyInputs;
use io::inputs::SimulationInputs;
use common::random::get_std_rng_for_rep;
use std::path::Path;
pub enum ChronologyExperiment{

    RandomWalk {experiment: RandomWalkSimulation },

    RJMCMC {experiment: RJMCMCSimulation },

}

impl ChronologyExperiment {

    pub fn setup( profile_inputs: &ChronologyInputs, 
                  inputs: &SimulationInputs, 
                  repetitions: usize, 
                  experiments: usize, 
                  minimum_traps:usize) 
        -> Result<Self,String>
    {

        if profile_inputs.setup.random_paths {
            let experiment = RandomWalkSimulation::new(profile_inputs.setup.fill_ratio.clone(), 
                                                                            profile_inputs.setup.fill_ratio_error.clone(), 
                                                                            profile_inputs.random_walks.clone(), 
                                                                            inputs, 
                                                                            repetitions,
                                                                            experiments, 
                                                                            minimum_traps)?;
            return Ok(Self::RandomWalk { experiment });
        }else{
            let experiment =  RJMCMCSimulation::new(profile_inputs.setup.fill_ratio.clone(), 
                                                                      profile_inputs.setup.fill_ratio_error.clone(), 
                                                                      profile_inputs.rjmcmc.clone(), 
                                                                      inputs, repetitions, 
                                                                      experiments, 
                                                                      minimum_traps)?;
            return Ok(Self::RJMCMC {experiment});
        }

    }

    pub fn run(&self, trial_num: usize, seed: usize) -> Result<(), String> {
        let mut rng = get_std_rng_for_rep(seed);
        match self {
            ChronologyExperiment::RandomWalk { experiment } => {
                experiment.run_random_walks(trial_num, Path::new("tmp"), & mut rng)?;
                Ok(())
            },
            ChronologyExperiment::RJMCMC { experiment } => {
                experiment.run_rjmcmc(trial_num, Path::new("tmp"), & mut rng)?;
                Ok(())
            }
        }

    }

}