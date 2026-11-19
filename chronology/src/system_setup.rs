// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Connect profile proposals, forward simulations, acceptance, and result storage.
//!
//! A run may tune a Bayesian walker first, then evaluate production proposals
//! and flush accepted profiles to the configured compressed output file.

use rand::Rng;
use std::path::Path;

use crate::profiles::{ProfileType, TimeTempProfile,MoveKind};
use crate::acceptance::{AcceptanceConstructor, AcceptorType};
use crate::burn_in::BurnInApparatus;
use common::constants::{temperature::TemperatureUnit,time::TimeUnit};
use common::random::get_std_rng_for_rep;
use io::chrono_inputs::ChronologyInputs;
use io::inputs::ResolvedExperiment;
use mc::system_setup::MonteCarloSimulation;


/// Profile generator, acceptance state, and forward Monte Carlo model for inference.
///
/// Coordinates are interpreted using the chronology time and temperature units.
pub struct ChronologyExperiment{
    generator:ProfileType ,
    acceptor:AcceptanceConstructor,
    monte_carlo:MonteCarloSimulation,
    time_unit:TimeUnit,
    temp_unit:TemperatureUnit,
}


impl ChronologyExperiment {

    /// Build the proposal generator, acceptor, and forward model from resolved inputs.
    ///
    /// `repetitions` controls independent trajectories per experiment; `experiments`
    /// selects configurations. Configuration and model construction errors propagate.
    pub fn setup( profile_inputs: &ChronologyInputs, 
                  inputs: &Vec<ResolvedExperiment>, 
                  repetitions: usize, 
                  experiments: usize, 
    ) -> Result<Self,String>
    {
        
        let generator = ProfileType::new(profile_inputs)?;
        let acceptor= AcceptanceConstructor::new(profile_inputs)?;
        let monte_carlo = MonteCarloSimulation::new(inputs, repetitions, experiments)?;

        
        Ok(Self{
            generator,
            acceptor,
            monte_carlo,
            time_unit: profile_inputs.setup.time_unit,
            temp_unit: profile_inputs.setup.temp_unit,
        })
        

    }
    
    /// Evaluate up to `trial_num` proposals and persist accepted profiles.
    ///
    /// Creates/replaces the chronology file in the existing output directory.
    /// Invalid proposals count toward the trial limit but skip forward simulation.
    /// The current Bayesian walk profile is updated by acceptance; remaining
    /// records are flushed on success. Generation, simulation, and I/O errors propagate.
    pub fn run_experiment(
        &mut self, 
        trial_num: usize, 
        output_directory: impl AsRef<Path>, 
        current_profile: &mut Option<TimeTempProfile>, 
        rng: &mut impl Rng
    ) -> Result<(), String> 
    {
        self.acceptor.set_output_file(output_directory, self.time_unit, self.temp_unit)?;

        let mut acct_cnt:usize = 0;

        for num in 0.. trial_num {
            
            let (trial_profile, move_kind, k_choice) = self.generator.generate_next_profile(&current_profile, rng)?;
            match trial_profile {
                None => {eprintln!("Profile {}/{trial_num} not acceptable | {acct_cnt} accepted! ",num+1,);
                        continue
                        },
                Some(trial_profile) =>{
                    let profile = trial_profile.generate_profile(self.time_unit)?;
                    let result = self.monte_carlo.run_to_final_ratio_only(profile, num)?;
                    acct_cnt += self.acceptor.accept_reject_profile(
                                                trial_profile,
                                                current_profile,
                                                move_kind,
                                                k_choice,
                                                result,
                                                self.monte_carlo.repetitions,
                                                rng)? as usize;
                    eprintln!("Checked profile number {}/{trial_num} | {acct_cnt} accepted! ",num+1,);

                }
            }
            
           
        }
        self.acceptor.final_write_check()?;

        Ok(())

    }

    /// Tune Bayesian walk proposal scales and return the final current profile.
    ///
    /// Requires a Bayesian walk setup and positive window/stability limits.
    /// Evaluates a starting profile, then adapts scales until convergence or the
    /// window limit. Copies tuned scales into the generator and acceptance rule.
    /// Invalid settings, failed initialization, and simulation errors are returned.
    /// No production output file is configured by this method.
    pub fn run_burn_in(
        &mut self, 
        inputs: &ChronologyInputs, 
        rng: &mut impl Rng
    ) -> Result<Option<TimeTempProfile>, String> 
    {
        let burn_in = inputs.burn_in;
        if burn_in.window == 0 || burn_in.maximum_no_windows == 0 || burn_in.stable_window == 0 {
            return Err("Burn-in window, maximum number of windows, and stable window count must all be greater than zero".to_string());
        }

        let mut burn = BurnInApparatus::new(inputs)?; 
       
        let mut current_profile:Option<TimeTempProfile> = None; 
    
        let (trial_profile, 
            _, 
            k_choice) = self.generator.generate_next_profile(&None, rng)?;
        
        match trial_profile {
            None =>  return Err("Burn-in failed unable to generate a starting profile".to_string()),
            Some(trial_profile) => {
                let profile= trial_profile.generate_profile(self.time_unit)?;
                let result = self.monte_carlo.run_to_final_ratio_only(profile, 0)?;
                let _ = self.acceptor.accept_reject_profile(
                        trial_profile,
                        &mut current_profile,
                        Some(MoveKind::Temp),
                        k_choice,
                        result,
                        self.monte_carlo.repetitions,
                        rng)?;

                    }
        }

        let mut mc_trial_num:usize = 1;
        for window_index in 0..burn_in.maximum_no_windows {
            for _ in 0..burn_in.window {
                let (trial_profile, move_kind, k_choice) = match &current_profile {
                    Some(current_profile) => {
                        burn.walker_alt.construct_new_profile(&current_profile, rng)?},
                    None => {return Err("There should now be a profile to use".to_string())}
                };
                let trial_move = burn.get_window(move_kind)?;
                if !self.generator.check_profile(&trial_profile){
                    trial_move.record(false);
                    continue
                }
                
                let profile = trial_profile.generate_profile(self.time_unit)?;
                let result = self.monte_carlo.run_to_final_ratio_only(profile, mc_trial_num)?;
             

                trial_move.record(self.acceptor.accept_reject_profile(
                                        trial_profile,
                                        &mut current_profile,
                                        Some(move_kind),
                                        k_choice,
                                        result,
                                        self.monte_carlo.repetitions, 
                                        rng)?);
               
                mc_trial_num += 1;
            }
            
            let converged = burn.burn_in_operation(window_index, 
                burn_in.maximum_no_windows, 
                burn_in.stable_window,)?;

            match &mut self.acceptor.acceptor {
                AcceptorType::BayesianWalk{values} => 
                    values.update_birth_sigmas(
                        burn.walker_alt.move_sigmas.time_birth*burn.walker_alt.move_sigmas.temp_birth
                    ),
                _ => unreachable!("Acceptor type should always be BayesianWalk for Burn in")
            };

            if converged {
                break;
            }
            burn.print_update(window_index + 1,burn_in.maximum_no_windows);
            
        }
        
       
        match &mut self.generator {
            ProfileType::Walk { walker, ..} => {
               walker.update_self(burn.walker_alt);
            },
            _ => unreachable!("Acceptor type should always be BayesianWalk for Burn in")
        };

        Ok(current_profile)
    }


    /// Seed profile generation, optionally tune a Bayesian walk, then run inference.
    ///
    /// `seed` identifies a stream derived from the global base seed. Burn-in
    /// runs only when enabled for Bayesian walk generation. Production output
    /// is written under the existing `tmp/` directory; errors propagate.
    pub fn run(
        &mut self, 
        inputs: &ChronologyInputs, 
        trial_num: usize, 
        seed: usize,
    ) -> Result<(), String> 
    {
        let mut rng = get_std_rng_for_rep(seed);
        
        let mut current_profile: Option<TimeTempProfile> = if 
        inputs.burn_in.run && inputs.setup.profile_generation == 'W' && inputs.setup.acceptance == 'B'{
            self.run_burn_in(&inputs, &mut rng)?

        }else { None};

        self.run_experiment(trial_num, Path::new("tmp"), &mut current_profile, &mut rng)
        
    }



}
