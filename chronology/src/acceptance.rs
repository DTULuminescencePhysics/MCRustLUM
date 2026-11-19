// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Compare simulated fills and ages with observations and record accepted histories.
//!
//! Residuals are normalized by measurement uncertainties. Bayesian rules use
//! log likelihoods; frequentist rules combine p-values. Accepted histories can
//! be buffered into compressed batches independently of the acceptance calculation.
use rand::Rng;
use statrs::distribution::{ChiSquared, ContinuousCDF};
use std::path::{Path, PathBuf};

use crate::profiles::{TimeTempProfile, MoveKind};
use common::numeric::{Float, TimeFloat};
use common::constants::{temperature::TemperatureUnit, time::TimeUnit};
use common::random::generate_log_open_uniform_random;
use io::chrono_inputs::ChronologyInputs;
use io::chrono_outputs::{
    AcceptanceInfo, AcceptedProfile,
    append_chronology_profile_to_experiment_file, 
    create_temporary_chronology_experiment_file, 
    write_temporary_chronology_units,};



/// Likelihood state and proposal scales used by the Bayesian profile walk.
///
/// Acceptance compares log likelihoods and includes a move-specific prior ratio.
/// The first proposal initializes the current profile without a likelihood comparison.
pub struct BayesianWalk{
    /// Gaussian log likelihood of the accepted/trial profile, without normalization constants.
    pub log_likelihood: Float,
    /// Absolute width of the allowed temperature prior range.
    pub temp_diff: Float,
    /// Product of the time and temperature birth scales used in prior corrections.
    pub birth_sigmas: Float, 

}

impl BayesianWalk {

    /// Store the log likelihood of the newly accepted current profile.
    pub fn update_likelihood(&mut self, new_likelihood: Float) {
        self.log_likelihood = new_likelihood;
    }
    /// Store the product of the time and temperature birth scales after tuning.
    pub fn update_birth_sigmas(&mut self, new_sigmas: Float) {
        self.birth_sigmas = new_sigmas;
    }

    /// Convert squared, uncertainty-normalised fill and age residuals to a Gaussian log likelihood.
    ///
    /// Returns `-0.5 * (sum(ratio_likelihood) + age_likelihood)`; normalization
    /// constants that cancel in likelihood ratios are omitted.
    fn calculate_log_likelihood(&self, ratio_likelihood: Vec<Float>, age_likelihood: Float) -> Float {
        -0.5*(ratio_likelihood.iter().sum::<Float>() + age_likelihood)
    }
    
    /// Calculate the log prior/proposal correction for inserting node `k_choice`.
    ///
    /// Uses the two new time intervals, the old interval, and the birth-scale product.
    /// The index must identify an interior node of the trial profile.
    fn birth_prior_ratio(
        &self,
        trial_profile: &TimeTempProfile, 
        current_profile: &TimeTempProfile, 
        k_choice: usize
    ) -> Float{
        (   
            (( trial_profile.times[k_choice-1]-trial_profile.times[k_choice])
                *(trial_profile.times[k_choice]-trial_profile.times[k_choice+1])
                * self.birth_sigmas)
            /((current_profile.times[k_choice-1]-current_profile.times[k_choice])
                *self.temp_diff)
            
        ).abs().ln()
    }

    /// Calculate the reverse birth correction for deleting an interior node.
    ///
    /// Endpoint deletions and profiles with no removable interior node return
    /// negative infinity so the proposal is rejected before indexing neighbours.
    fn death_prior_ratio(  
        &self,
        trial_profile: &TimeTempProfile, 
        current_profile: &TimeTempProfile, 
        k_choice: usize
    ) -> Float{
        // Birth only splits an existing interval, so it cannot reverse an
        // endpoint deletion. Reject these proposals before indexing neighbours.
        if current_profile.times.len() <= 2 || k_choice == 0
            || k_choice >= current_profile.times.len() - 1
        {
            return Float::NEG_INFINITY;
        }
        (
            ( (trial_profile.times[k_choice-1] - trial_profile.times[k_choice]) 
             * self.temp_diff)
            /( (current_profile.times[k_choice-1]-current_profile.times[k_choice])
              *(current_profile.times[k_choice]-current_profile.times[k_choice+1])
              *self.birth_sigmas)
        ).abs().ln()
    }

    /// Calculate the log prior ratio from the intervals adjacent to a moved time node.
    ///
    /// Endpoint moves use one interval; interior moves use the product of both.
    fn time_prior_ratio(
        &self,
        trial_profile: &TimeTempProfile, 
        current_profile: &TimeTempProfile, 
        k_choice: usize
    ) -> Float {
        if k_choice == 0 {
            (
                    (trial_profile.times[k_choice+1] - trial_profile.times[k_choice])
                /(current_profile.times[k_choice+1] - current_profile.times[k_choice])
            ).abs().ln()

        } else if k_choice == trial_profile.times.len() -1 {
            (
                    (trial_profile.times[k_choice] - trial_profile.times[k_choice-1])
                /(current_profile.times[k_choice] - current_profile.times[k_choice-1])
            ).abs().ln()

        }else{
            (
                    ( (trial_profile.times[k_choice] - trial_profile.times[k_choice-1])
                    *(trial_profile.times[k_choice+1] - trial_profile.times[k_choice]))
                /( (current_profile.times[k_choice] - current_profile.times[k_choice-1])
                    *(current_profile.times[k_choice+1] - current_profile.times[k_choice]))
            ).abs().ln()
        
        }
    }

    /// Apply a Bayesian walk decision and return acceptance plus diagnostic metadata.
    ///
    /// On acceptance, update the stored likelihood and the affected current-profile
    /// coordinates. A missing current profile is initialized from the trial.
    /// `k_choice` identifies the inserted, deleted, or perturbed node.
    /// The stored posterior field is currently a zero placeholder.
    fn calculate_acceptance_criteria(
        &mut self, 
        ratio_likelihood: Vec<Float>, 
        age_likelihood: Float,
        trial_profile: &TimeTempProfile, 
        current_profile: &mut Option<TimeTempProfile>,
        move_kind: MoveKind,
        k_choice: usize, 
        rng: &mut impl Rng
    )-> (bool, AcceptanceInfo) {

        let log_likelihood = self.calculate_log_likelihood(ratio_likelihood, age_likelihood);
        match current_profile {
            None => {
                self.update_likelihood(log_likelihood);
                *current_profile = Some(trial_profile.clone());
                
                return (true, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 });
            }, 
            Some(current_profile) => {
               
                match move_kind {
                    MoveKind::Birth => {
                        let prior_ratio = self.birth_prior_ratio(trial_profile, current_profile, k_choice);
                        if generate_log_open_uniform_random(rng) < log_likelihood - self.log_likelihood + prior_ratio{
                            self.update_likelihood(log_likelihood);
                            *current_profile = trial_profile.clone();
                            (true, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        } else {
                            (false, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        }

                    },
                    MoveKind::Death => {
                        let prior_ratio = self.death_prior_ratio(trial_profile, current_profile, k_choice);
                        if generate_log_open_uniform_random(rng) < log_likelihood - self.log_likelihood + prior_ratio{
                            self.update_likelihood(log_likelihood);
                            *current_profile = trial_profile.clone();
                            (true, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        } else {
                           (false, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        }

                    },
                    MoveKind::Time => {
                        let prior_ratio = self.time_prior_ratio(trial_profile, current_profile, k_choice);
                        if generate_log_open_uniform_random(rng) < log_likelihood - self.log_likelihood + prior_ratio{
                            self.update_likelihood(log_likelihood);
                            current_profile.times = trial_profile.times.clone();
                            (true, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        } else {
                           (false, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        }
                    },
                    MoveKind::Temp => {
                        if generate_log_open_uniform_random(rng) < log_likelihood - self.log_likelihood{
                            self.update_likelihood(log_likelihood);
                            current_profile.temps = trial_profile.temps.clone();
                            (true, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        } else {
                           (false, AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior: 0.0 })
                        }
                    }
                }
                   
            }

        }

    }

}



/// Acceptance rule paired with independent random profiles or a profile walk.
///
/// Frequentist modes use a combined p-value with a 0.05 rejection threshold;
/// Bayesian modes compare a uniform log draw with a log-likelihood criterion.
pub enum AcceptorType{
    /// Independent rejection sampling using the trial log likelihood.
    BayesianRandom,
    /// Bayesian walk acceptance with move-specific prior/proposal corrections.
    BayesianWalk{
        /// Current Bayesian likelihood and proposal-correction state.
        values: BayesianWalk,
    },
    /// Independent profile acceptance using a combined frequentist p-value.
    FrequentistRandom,
    /// Walk proposal testing using a combined frequentist p-value.
    FrequentistWalk,
}

impl AcceptorType{

    /// Select an acceptance strategy from the chronology setup flags.
    ///
    /// Returns an error unless acceptance is `F` or `B` and generation is `R` or `W`.
    pub fn new(inputs: &ChronologyInputs) -> Result<Self, String>{
       
        match inputs.setup.acceptance{
            'F' => {
                match inputs.setup.profile_generation {
                    'R' => Ok(AcceptorType::FrequentistRandom),
                    'W' => Ok(AcceptorType::FrequentistWalk),
                    _   => Err(format!("Invalid profile generation character: only R, and W are allowed, {} was entered",inputs.setup.profile_generation))
                }
            },
            'B' => {
                match inputs.setup.profile_generation {
                    'R' => Ok(AcceptorType::BayesianRandom),
                    'W' => Ok(AcceptorType::BayesianWalk{ 
                        values: BayesianWalk { 
                            log_likelihood: Float::INFINITY, 
                            temp_diff: (inputs.constraints.temp_range.0-inputs.constraints.temp_range.1).abs(), 
                            birth_sigmas: inputs.random_walk.s_temp_birth*inputs.random_walk.s_time_birth
                         }}),
                    _   => Err(format!("Invalid profile generation character: only R, and W are allowed, {} was entered",inputs.setup.profile_generation))
                }

            },
            _ => {
                Err("Invalid profile acceptance character: only F, and B are allowed".into())
            }
        }
    }
    /// Return the upper-tail chi-squared probability for a squared residual.
    ///
    /// Rejects zero degrees of freedom and negative or non-finite residuals.
    fn chi_squared_test(&self,
        squared_residual: &Float,
        degrees_of_freedom: usize,
    ) -> Result<Float, String> 
    {
        if degrees_of_freedom == 0 {
            return Err("Degrees of freedom must be positive".into());
        }
         if !squared_residual.is_finite()  {
            return Err(format!("Squared residuals must be finite, but had value {squared_residual}"));
        }
        if *squared_residual < 0.0 {
            return Err(format!("Squared residuals must be nonnegative, but had value {squared_residual}"))
        }

        let distribution = ChiSquared::new(degrees_of_freedom as f64)
            .map_err(|error| error.to_string())?;

        let p_value = distribution.sf(*squared_residual);

        Ok(p_value)
    }

    /// Combine individual fill-ratio and age p-values using Fisher's method.
    ///
    /// Each squared residual is tested with one degree of freedom; the combined
    /// statistic is tested with twice the number of measurements. Invalid
    /// residuals propagate as errors.
    pub fn frequentist_testing(
        &self,  
        ratio_likelihood: Vec<Float>, 
        age_likelihood: Float,
    ) -> Result<Float, String> 
    {
        let n = ratio_likelihood.len()+1;
        let mut p_values: Vec<Float> = Vec::with_capacity(n);
        for ll in ratio_likelihood.iter(){
            p_values.push(self.chi_squared_test(ll,1)?);
        }
        p_values.push(self.chi_squared_test(&age_likelihood,1)?);
        let summed: Float = -2.0 * p_values.iter().map(|p| p.ln()).sum::<Float>();
        
        self.chi_squared_test(&summed, 2*n)
    
    }

    /// Evaluate a trial and return `(accepted, metadata)` for the selected strategy.
    ///
    /// Only the Bayesian walk updates `current_profile` here. That mode rejects
    /// a trial without a move kind. Statistical validation errors are returned.
    pub fn make_decision(
        &mut self, 
        ratio_likelihood: Vec<Float>, 
        age_likelihood: Float,
        trial_profile: &TimeTempProfile,
        current_profile: &mut Option<TimeTempProfile>,
        move_kind: Option<MoveKind>, 
        k_choice: usize,
        rng: &mut impl Rng
    ) -> Result<(bool, AcceptanceInfo), String> 
    {   
        match self {
            AcceptorType::BayesianRandom => {
                let log_likelihood = -0.5*(ratio_likelihood.iter().sum::<Float>() + age_likelihood);
                if generate_log_open_uniform_random(rng) < log_likelihood{
                    return Ok((true, AcceptanceInfo::BayesianRandom { log_likelihood })); 
                }else{
                    return Ok((false, AcceptanceInfo::BayesianRandom { log_likelihood }));
                }
            },
            AcceptorType::BayesianWalk {values} => {
                match move_kind {
                    Some(move_kind) => {
                        let result = values.calculate_acceptance_criteria(
                            ratio_likelihood, 
                            age_likelihood, 
                            trial_profile, 
                            current_profile, 
                            move_kind, 
                            k_choice, 
                            rng); 
                        return Ok(result)
                    },
                    None => return Ok((false, AcceptanceInfo::BayesianWalk { log_likelihood: 0.0, log_posterior: 0.0 }))
                }

            },
            AcceptorType::FrequentistRandom => {
                let full_p = self.frequentist_testing(ratio_likelihood,age_likelihood)?;
               
                if full_p < 0.05 {
                    return Ok((false, AcceptanceInfo::Frequentist { p_test: full_p })); 
                }else {
                    return Ok((true, AcceptanceInfo::Frequentist { p_test: full_p })); 
                }
            },
            AcceptorType::FrequentistWalk => {
                let full_p = self.frequentist_testing(ratio_likelihood,age_likelihood)?;
               
                if full_p < 0.05 {
                    return Ok((false, AcceptanceInfo::Frequentist { p_test: full_p })); 
                }else {
                    return Ok((true, AcceptanceInfo::Frequentist { p_test: full_p })); 
                }

            }

        }
    }

}


/// Measured targets, acceptance state, and a buffer of accepted profiles.
///
/// Fill measurements correspond to experiment order. Accepted trials are
/// written in batches of 1,000 when an output file has been configured.
pub struct AcceptanceConstructor {
    /// Measured occupied-trap fraction for each experiment, in experiment order.
    pub fill_ratio: Vec<Float>,
    /// Non-zero measurement uncertainty for each corresponding fill ratio.
    pub fill_ratio_error: Vec<Float>,
    /// Measured profile duration in the setup time unit.
    pub age: TimeFloat, 
    /// Non-zero age uncertainty in the setup time unit.
    pub age_error:TimeFloat,
    /// Acceptance strategy and its current statistical state.
    pub acceptor: AcceptorType,
    /// Optional output buffer of accepted trials; `None` disables persistence.
    pub accepted_profiles: Option<Vec<AcceptedProfile<TimeTempProfile>>>,
    /// Destination of the compressed accepted-profile batches.
    pub output_path: PathBuf, 

}

impl AcceptanceConstructor{

    /// Copy measured targets and uncertainties and construct the requested acceptor.
    ///
    /// The age and age-error arrays must each contain at least one value; only
    /// the first age measurement is used. Output remains disabled until
    /// `set_output_file` is called. Invalid acceptance flags return an error.
    pub fn new(inputs: &ChronologyInputs) -> Result<Self, String>{
        Ok(Self{
            fill_ratio:inputs.measurement.fill_ratio.clone(),
            fill_ratio_error:inputs.measurement.fill_ratio_error.clone(),
            age: inputs.measurement.age[0],
            age_error: inputs.measurement.age_error[0],
            acceptor: AcceptorType::new(inputs)?,
            accepted_profiles: None,
            output_path:  PathBuf::default(),
        })
    }

    /// Return squared normalized residuals between mean simulated and measured fills.
    ///
    /// Each inner result vector contains repetitions for one experiment. `reps`
    /// and measurement uncertainties must be non-zero. Zipping limits the
    /// comparison to the shortest of the result, target, and uncertainty vectors.
    fn ratio_likelihood(&self, result: Vec<Vec<Float>>, reps: usize) -> Vec<Float> {

        result.iter()
              .zip(&self.fill_ratio)
              .zip(&self.fill_ratio_error)
              .map(|((res, target), error)|{
                   let diff: Float = (((res.iter().sum::<Float>())/reps as Float)-target)/error;
                   diff.powi(2)}
              )
              .collect() 
    }

    /// Return the squared age residual divided by the squared age uncertainty.
    ///
    /// The trial age and measured age must use the same unit and the uncertainty
    /// must be non-zero.
    fn age_likelihood(&self, age:TimeFloat) -> Float {
        ((age-self.age)/self.age_error).powi(2) as Float
    }

    /// Score a simulated trial, apply the selected rule, and buffer accepted trials.
    ///
    /// `result` is indexed by experiment then repetition. Returns whether the trial
    /// was accepted, or a statistical/output error. Without an output buffer
    /// (e.g. during burn-in), acceptance is evaluated without recording the trial.
    pub fn accept_reject_profile(&mut self, 
        trial_profile: TimeTempProfile, 
        current_profile: &mut Option<TimeTempProfile>,
        move_kind: Option<MoveKind>,
        k_choice: usize,
        result: Vec<Vec<Float>>, 
        reps: usize,
        rng: &mut impl Rng
    ) -> Result<bool,String> 
    {

        let ratio_likelihood = self.ratio_likelihood(result, reps);
        let age_likelihood = self.age_likelihood(trial_profile.age());
        let (add_prof, data) =self.acceptor.make_decision(
            ratio_likelihood, 
            age_likelihood, 
            &trial_profile, 
            current_profile, 
            move_kind, 
            k_choice, 
            rng)?; 

        if add_prof {
            if let Some(accepted_profiles) = self.accepted_profiles.as_mut() {
                accepted_profiles.push(AcceptedProfile { profile: trial_profile, acceptance: data });
                if accepted_profiles.len() >= 1000 {
                    self.write_to_output()?;
                }
            }
        }
        Ok(add_prof)
    }
    /// Flush any remaining accepted profiles at the end of a run.
    pub fn final_write_check(&mut self) -> Result<(), String> {
        self.write_to_output()
    }

    /// Append the buffered accepted profiles and clear them only after a successful write.
    ///
    /// An absent or empty buffer is a no-op. A failed write returns an error and
    /// retains the buffer for inspection or retry.
    pub fn write_to_output(&mut self) -> Result<(), String> {
        if let Some(accepted_profiles) = self.accepted_profiles.as_mut() {
            if !accepted_profiles.is_empty() {
                append_chronology_profile_to_experiment_file(
                    &self.output_path, accepted_profiles)
                    .map_err(|error| error.to_string())?;
                accepted_profiles.clear();
            }
        }
        Ok(())

    }

    /// Create `chronology_results.bin.gz`, write its unit header, and enable buffering.
    ///
    /// The directory must already exist. Existing contents at that path are
    /// replaced. The units describe the profile coordinates stored in later batches;
    /// file creation and header-writing errors are returned.
    pub fn set_output_file(
        &mut self, 
        output_directory: impl AsRef<Path>, 
        time_unit: TimeUnit, 
        temp_unit: TemperatureUnit,
    )-> Result<(), String> 
    {
        let output_directory = output_directory.as_ref();
        self.output_path = output_directory.join("chronology_results.bin.gz");
        
        create_temporary_chronology_experiment_file(
            &self.output_path)
            .map_err(|error| error.to_string())?;
           
        write_temporary_chronology_units(
            &self.output_path, 
            time_unit, 
            temp_unit)
            .map_err(|error| error.to_string())?;

        self.accepted_profiles = Some(Vec::with_capacity(1000));

        Ok(())

    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use io::chrono_outputs::read_all_chrono_batches;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    /// Verify that writes full batches and flushes only remaining accepted profiles.
    fn writes_full_batches_and_flushes_only_remaining_accepted_profiles() {
        let unique = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "mcrustlum_acceptance_{}_{unique}", std::process::id()));
        fs::create_dir(&directory).unwrap();
        let mut inputs = ChronologyInputs::default();
        inputs.measurement.age = vec![1.0];
        let mut acceptor = AcceptanceConstructor::new(&inputs).unwrap();
        let mut rng = common::random::get_std_rng_for_rep(0);
        let mut current = None;
        let profile = TimeTempProfile {
            times: vec![0.0, 1.0],
            temps: vec![100.0, 110.0],
        };

        // Burn-in can make decisions before an output file is configured.
        assert!(acceptor.accept_reject_profile(profile.clone(), &mut current,
            None, 0, vec![vec![0.0]], 1, &mut rng).unwrap());
        acceptor.final_write_check().unwrap();
        assert!(acceptor.accepted_profiles.is_none());

        acceptor.set_output_file(&directory, TimeUnit::Second, TemperatureUnit::Celsius).unwrap();
        // A rejected trial must not enter the output buffer.
        assert!(!acceptor.accept_reject_profile(profile.clone(), &mut current,
            None, 0, vec![vec![1.0]], 1, &mut rng).unwrap());
        for _ in 0..1001 {
            assert!(acceptor.accept_reject_profile(profile.clone(), &mut current,
                None, 0, vec![vec![0.0]], 1, &mut rng).unwrap());
        }
        assert_eq!(acceptor.accepted_profiles.as_ref().unwrap().len(), 1);
        acceptor.final_write_check().unwrap();
        acceptor.final_write_check().unwrap();
        let reader = read_all_chrono_batches::<AcceptedProfile<TimeTempProfile>>(&acceptor.output_path).unwrap();
        assert_eq!(reader.time_unit, TimeUnit::Second);
        assert_eq!(reader.temp_unit, TemperatureUnit::Celsius);
        let batches = reader.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(batches.iter().map(Vec::len).collect::<Vec<_>>(), vec![1000, 1]);
        for stored in batches.iter().flatten() {
            assert_eq!(stored.profile.times, profile.times);
            assert_eq!(stored.profile.temps, profile.temps);
            assert!(matches!(stored.acceptance,
                AcceptanceInfo::Frequentist { p_test: 1.0 }));
        }
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    /// Verify that failed write keeps buffered profiles.
    fn failed_write_keeps_buffered_profiles() {
        let mut acceptor = AcceptanceConstructor::new(&ChronologyInputs::default()).unwrap();
        acceptor.accepted_profiles = Some(vec![ AcceptedProfile{profile: TimeTempProfile {
            times: vec![0.0, 1.0], temps: vec![100.0, 110.0]}, acceptance: AcceptanceInfo::BayesianRandom {log_likelihood: 0.5},
        }]);
        // The default empty output path cannot be opened as a file.
        assert!(acceptor.final_write_check().is_err());
        assert_eq!(acceptor.accepted_profiles.as_ref().unwrap().len(), 1);
    }

    #[test]
    /// Verify that stores acceptance metadata for every acceptor type.
    fn stores_acceptance_metadata_for_every_acceptor_type() {
        let unique = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "mcrustlum_metadata_{}_{unique}", std::process::id()));
        fs::create_dir(&directory).unwrap();
        let profile = TimeTempProfile {
            times: vec![0.0, 1.0],
            temps: vec![100.0, 110.0],
        };

        for (acceptance, generation) in [('B', 'R'), ('B', 'W'), ('F', 'R'), ('F', 'W')] {
            let mut inputs = ChronologyInputs::default();
            inputs.setup.acceptance = acceptance;
            inputs.setup.profile_generation = generation;
            inputs.measurement.age = vec![1.0];
            inputs.measurement.age_error = vec![1.0];
            inputs.measurement.fill_ratio = vec![0.0];
            inputs.measurement.fill_ratio_error = vec![1.0];
            let mut acceptor = AcceptanceConstructor::new(&inputs).unwrap();
            acceptor.set_output_file(&directory, TimeUnit::Second, TemperatureUnit::Celsius).unwrap();
            let mut current = None;
            let mut rng = common::random::get_std_rng_for_rep(0);
            // The residual is 0.5^2, giving a nonzero log likelihood of -0.125.
            let accepted = (0..100).any(|_| acceptor.accept_reject_profile(
                profile.clone(), &mut current, Some(MoveKind::Temp), 0,
                vec![vec![0.5]], 1, &mut rng,
            ).unwrap());
            assert!(accepted);
            acceptor.final_write_check().unwrap();
            let batches = read_all_chrono_batches::<AcceptedProfile<TimeTempProfile>>(
                &acceptor.output_path,
            ).unwrap().collect::<Result<Vec<_>, _>>().unwrap();
            assert_eq!(batches.len(), 1);
            assert_eq!(batches[0].len(), 1);
            let stored = &batches[0][0];
            assert_eq!(stored.profile.times, profile.times);
            assert_eq!(stored.profile.temps, profile.temps);
            match (acceptance, generation, &stored.acceptance) {
                ('B', 'R', AcceptanceInfo::BayesianRandom { log_likelihood }) => {
                    assert_eq!(*log_likelihood, -0.125);
                }
                ('B', 'W', AcceptanceInfo::BayesianWalk { log_likelihood, log_posterior }) => {
                    assert_eq!(*log_likelihood, -0.125);
                    assert_eq!(*log_posterior, 0.0);
                }
                ('F', _, AcceptanceInfo::Frequentist { p_test }) => {
                    // chi-square(1) survival at 0.25, combined with an age p-value of 1.
                    let p: Float = 0.6170750774519738;
                    let expected = p * (1.0 - p.ln());
                    assert!((*p_test - expected).abs() < 1.0e-12);
                }
                _ => panic!("unexpected metadata for {acceptance}/{generation}"),
            }
        }
        fs::remove_dir_all(directory).unwrap();
    }
}
