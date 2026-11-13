// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Creation and per-iteration reset of Monte Carlo state.

use crate::experiment::MCExperiment;
use common::crystal::Cube;
use common::numeric::Float;
use common::rate_equation_selection::Transitions;
use common::time_temperature::TimeTemperature;
use io::inputs::{CubeDefinition, 
    ResolvedExperiment, TimeTemperatureDefinition, TrapDefinition, BandDefinition, 
    FillingDefinition, TransitionDefinition};
use common::trap_hole_band_tail::{TrapGenerator, ReTrapParameterLayout};
use common::rate_equation_inputs::FillingTransitionInputs;
use std::path::Path;
use rayon::prelude::*;
// use std::sync::atomic::{AtomicUsize, Ordering};

/// State shared by the repetitions and experiments in one Monte Carlo run.
///
/// Constructing this type validates the supplied geometry and temperature
/// profile and randomly places the crystal sites.
#[derive(Debug, Clone)]
pub struct MonteCarloSimulation {
    /// Actual coordinates are generated independently for each repetition.
    pub cubes: Vec<Cube>,
    /// Current time and temperature state.
    pub time_temperature: Vec<TimeTemperature>,
    /// Enabled physical pathways and their rate equations.
    pub transitions: Vec<Transitions>,
    /// Vector containing trap generators that are used to 
    /// set the physical properties of the electron traps for each experiment
    pub trap_generators: Vec<Vec<TrapGenerator>>,
    /// Vector containing the filling information 
    pub filling: Vec<Option<FillingTransitionInputs>>,
    /// Vector holding retrapping parameters for each experiment
    pub retrapping: Vec<ReTrapParameterLayout>,
    /// Vector holding the initial number of occupied traps
    pub initial_traps: Vec<Vec<usize>>,
    /// Vector holding the initial number of occupied holes
    pub initial_holes:Vec<usize>,
    /// Number of independent stochastic repetitions performed per experiment.
    /// The field name preserves the crate's existing `repetitions` spelling.
    pub repetitions: usize,
    /// Number of independent experiments to run.
    pub experiments: usize,
}

impl MonteCarloSimulation {
    /// Build all runtime state from a grouped input configuration.
    
    pub fn new(
        inputs: &Vec<ResolvedExperiment>,
        repetitions: usize,
        experiments: usize,
        minimum_traps:usize,
    ) -> Result<Self, String> {
        if experiments == 0 {
            return Err("experiments must be positive".to_string());
        }
        
        let cubes = (0..experiments)
            .map(|experiment_index| {
                let densities: Vec<Float> = inputs[experiment_index]
                    .traps
                    .iter().map(|(_, trap)| trap.density.value as Float)
                    .collect();
                   
                MonteCarloSimulation::generate_cube_for_experiment(
                    &inputs[experiment_index].cube,
                    densities,
                    minimum_traps,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let initial_traps: Vec<Vec<usize>> = (0..experiments)
            .map(|experiment_index| {
                cubes[experiment_index].trap_total
                .iter()
                .zip(&inputs[experiment_index].traps)
                .map(|(tot,(_,trap))| {
                    (*tot as Float * trap.initial_fill.value).round() as usize
                }).collect()
            }).collect();
        
        let initial_holes: Vec<usize> = (0..experiments)
            .map(|experiment_index| 
                initial_traps[experiment_index].iter().sum()
            ).collect();

        let time_temperature = (0..experiments)
            .map(|experiment_index| {
                MonteCarloSimulation::generate_time_temperature(&inputs[experiment_index].time_temperature)
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        let transitions = (0..experiments)
            .map(|experiment_index| {
                MonteCarloSimulation::generate_transitions(&inputs[experiment_index].transitions)
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        let trap_generators = (0..experiments)
            .map(|experiment_index| {
                MonteCarloSimulation::generate_trap_generators( &inputs[experiment_index].traps)
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        let filling = (0..experiments)
            .map(|experiment_index| {
                let total_traps: usize = cubes[experiment_index].trap_total.iter().sum();
                MonteCarloSimulation::generate_filling(
                    &inputs[experiment_index].filling, 
                    inputs[experiment_index].transitions.filling,
                    total_traps)
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        let retrapping = (0..experiments)
        .map(|experiment_index| {
            MonteCarloSimulation::generate_retrapping(&inputs[experiment_index].transitions, &inputs[experiment_index].band)
        }).collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            cubes,
            time_temperature,
            transitions,
            trap_generators,
            filling,
            retrapping,
            initial_traps,
            initial_holes,
            repetitions,
            experiments,
        })
    }

    /// Generate the simulation boundary and derive its site counts.
    ///
    /// The trap count is derived from the configured density and volume. Hole
    /// and bandtail counts are interpreted as ratios per trap. If this would
    /// create fewer than `minimum_traps` traps, the volume is enlarged isotropically
    /// while preserving the requested density.
    /// Generate a cube using the density configured for one experiment.
    pub fn generate_cube_for_experiment(
        inputs: &CubeDefinition,
        densities: Vec<Float>,
        minimum_traps: usize,
    ) -> Result<Cube, String> {

        Cube::new_from_density(
            inputs.cube_height.value,
            inputs.cube_width.value,
            inputs.cube_length.value,
            densities,
            inputs.hole_count,
            inputs.bandtail_count,
            inputs.periodic,
            minimum_traps,
        )
    }
    /// Validate and generate the piecewise-linear time/temperature profile.
    pub fn generate_time_temperature(
        inputs: &TimeTemperatureDefinition,
    ) -> Result<TimeTemperature, String> {

        TimeTemperature::new(
            inputs.times.value.clone(),
            inputs.temperatures.value.clone(),
            inputs.times.unit.get_time_unit()?,
        )
    }
    /// Generate the transition selection used by the simulation.
    ///
    /// Delocalised release currently uses first-order kinetics. Boolean flags
    /// from the localised, delocalised, and filling input groups determine
    /// which pathways are active.
    pub fn generate_transitions(inputs: &TransitionDefinition) -> Result<Transitions, String> {
        Transitions::from_bool(
            inputs.ground_state_tunnelling,
            inputs.excited_state_tunnelling,
            inputs.ground_state_to_conduction_band,
            inputs.excited_state_to_conduction_band,
            inputs.filling,
            "first",
            inputs.ground_state_tunnelling_retrapping,
            inputs.excited_state_tunnelling_retrapping,
        )
    }
    /// Generate the filling selection used by the simulation
    pub fn generate_filling(inputs: &FillingDefinition, run_bool: bool, total_n: usize) -> Result<Option<FillingTransitionInputs>, String> {
        if run_bool{
            let fill = FillingTransitionInputs::get_inputs(
                       inputs.d0.value,
                    inputs.d_dot.value,
                        0,
                        total_n,
                    )?;
            Ok(Some(fill))
        
        } else {
            Ok(None)
        }
    }
    /// Generate the trap generators 
    pub fn generate_trap_generators( trap_information: &Vec<(String, TrapDefinition)> )-> Result<Vec<TrapGenerator>, String> {

       Ok(trap_information
            .iter()
            .map(|(name,trap_info)|{
                TrapGenerator::new(
                    name.clone(), 
                    trap_info.ground_to_excited_state_energy.value, 
                    trap_info.s_ground_to_excited_frequency.value, 
                    trap_info.s_excited_to_ground_frequency.value, 
                    trap_info.ground_to_conduction_band_energy.value,
                    trap_info.ground_to_conduction_band_energy.value-trap_info.ground_to_excited_state_energy.value, 
                    trap_info.s_ground_state_to_conduction_band_frequency.value, 
                    trap_info.s_excited_state_to_conduction_band_frequency.value, 
                    trap_info.b_ground_state_to_hole_frequency.value, 
                    trap_info.b_excited_state_to_hole_frequency.value, 
                    trap_info.alpha_ground_state_barrier_constant.value, 
                    trap_info.alpha_excited_state_barrier_constant.value, 
                    trap_info.ground_to_excited_state_energy.function.clone(), 
                    trap_info.s_ground_to_excited_frequency.function.clone(), 
                    trap_info.s_excited_to_ground_frequency.function.clone(), 
                    trap_info.ground_to_conduction_band_energy.function.clone(),
                    trap_info.s_ground_state_to_conduction_band_frequency.function.clone(), 
                    trap_info.s_excited_state_to_conduction_band_frequency.function.clone(), 
                    trap_info.b_ground_state_to_hole_frequency.function.clone(), 
                    trap_info.b_excited_state_to_hole_frequency.function.clone(),
                    trap_info.alpha_ground_state_barrier_constant.function.clone(), 
                    trap_info.alpha_excited_state_barrier_constant.function.clone(), 
                )
            })
            .collect()
        )
    }

    pub fn generate_retrapping(inputs: &TransitionDefinition, band: &BandDefinition) -> Result<ReTrapParameterLayout, String> {

        Ok(ReTrapParameterLayout::new(
            inputs.ground_state_to_conduction_band || inputs.excited_state_to_conduction_band, 
            inputs.filling, 
            inputs.conduction_band_retrapping, 
            inputs.filling_retrapping, 
            inputs.gaussian_kernel, 
            band.conduction_band_hole_to_trap_preference.value, 
            band.valence_band_hole_to_trap_preference.value, 
            band.mu_conduction_band_distance.value, 
            band.mu_valence_band_distance.value))

    }

    /// Restore the time/temperature profile to its first control point.
    pub fn reset_time_temperature(&mut self,i:usize) {
        self.time_temperature[i].reset();
    }

    /// Run every configured experiment and repetition.
    ///
    /// Every repetition receives a newly generated random spatial realization
    /// together with a reset copy of the time/temperature profile.
    pub fn run(&mut self) -> Result<(), String> {
        self.run_to_directory(Path::new("tmp"))
    }

    /// Run every configured experiment and repetition, placing temporary
    /// result files in `output_directory`.
    pub fn run_to_directory(&mut self, output_directory: impl AsRef<Path>) -> Result<(), String> {
        let output_directory = output_directory.as_ref();
        let batch_capacity: usize = 100;

        // Keep experiments sequential. Only the repetitions for the current
        // experiment are distributed across the Rayon thread pool.
        for experiment_index in 0..self.experiments {
            let cube = &self.cubes[experiment_index];
            

            (0..self.repetitions)
                .into_par_iter()
                .try_for_each(|repetition_index| {
                    let output_path = output_directory.join(format!(
                        "experiment_results_{experiment_index}_{repetition_index}.bin.gz"
                    ));

                    let mut time_temperature = self.time_temperature[experiment_index].clone();
                    time_temperature.reset();

                    // Unique across all experiment/repetition pairs.
                    let random_repetition = experiment_index * self.repetitions + repetition_index;
                    let mut filling = self.filling[experiment_index];

                    let mut experiment = MCExperiment::initialise(
                        cube,
                        &self.initial_traps[experiment_index],
                        &self.initial_holes[experiment_index],
                        &self.trap_generators[experiment_index],
                        self.retrapping[experiment_index],
                        &mut filling,
                        time_temperature,
                        &random_repetition,
                    )
                    .map_err(|error| {
                        format!(
                            "could not initialise experiment {} of {}, repetition {} of {}: {error}",
                            experiment_index + 1,
                            self.experiments,
                            repetition_index + 1,
                            self.repetitions,
                        )
                    })?;
                   

                    experiment
                        .run(
                            cube,
                            &self.transitions[experiment_index],
                            &output_path,
                            &batch_capacity,
                        )
                        .map_err(|error| {
                            format!(
                                "experiment {} of {}, repetition {} of {} failed: {error}",
                                experiment_index + 1,
                                self.experiments,
                                repetition_index + 1,
                                self.repetitions,
                            )
                        })
                })?;
        }

        Ok(())
    }

    /// Run every configured experiment and repetition and return final ratios.
    ///
    /// The outer vector is indexed by experiment. Each inner vector contains
    /// one final fill ratio per repetition; repetition order is unspecified.
    pub fn run_to_final_ratio_only(&self, tt_profile: TimeTemperature, trial_num:usize) -> Result<Vec<Vec<Float>>, String> {
        // let completed_repetitions = AtomicUsize::new(0);
        // let total_repetitions = self.experiments * self.repetitions;

        let mut results_by_experiment = Vec::with_capacity(self.experiments);
        for experiment_index in 0..self.experiments {
            let cube = &self.cubes[experiment_index];
           
            (0..self.repetitions)
                    .into_par_iter()
                    .map(|repetition_index| {

                        let mut time_temperature = tt_profile.clone();
                        time_temperature.reset();
                        let mut filling = self.filling[experiment_index];
                        // Unique across all experiment/repetition pairs.
                        let random_repetition =
                            (trial_num*(self.experiments*self.repetitions))+
                            experiment_index * self.repetitions + repetition_index;

                        let mut experiment = MCExperiment::initialise(
                            cube,
                            &self.initial_traps[experiment_index],
                            &self.initial_holes[experiment_index],
                            &self.trap_generators[experiment_index],
                            self.retrapping[experiment_index],
                            &mut filling,
                            time_temperature,
                            &random_repetition,
                        )
                        .map_err(|error| {
                        format!(
                            "could not initialise experiment {} of {}, repetition {} of {}: {error}",
                            experiment_index + 1,
                            self.experiments,
                            repetition_index + 1,
                            self.repetitions,
                        )
                        })?;
                        

                        let final_ratio = experiment
                            .run_no_output_file(
                                cube,
                                &self.transitions[experiment_index],
                            )
                            .map_err(|error| {
                                format!(
                                    "experiment {} of {}, repetition {} of {} failed: {error}",
                                    experiment_index + 1,
                                    self.experiments,
                                    repetition_index + 1,
                                    self.repetitions,
                                )
                            })?;

                        // let completed = completed_repetitions.fetch_add(1, Ordering::Relaxed) + 1;
                        // eprintln!("Repetitions completed: {completed}/{total_repetitions}");

                        Ok(final_ratio)
                    })
                    .collect::<Result<Vec<Float>, String>>()
                    .map(|experiment_results| results_by_experiment.push(experiment_results))?;
        }

        Ok(results_by_experiment)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::rate_equation_selection::{
        DelocalisedRateEquation, DelocalisedRateEquationType, FillingRateEquation,
        LocalisedRateEquation,
    };
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_output_directory(label: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_{label}_{}_{}",
            std::process::id(),
            unique,
        ))
    }

    fn small_inputs() -> Vec<ResolvedExperiment> {
        let fixture = include_str!("../../input.toml")
            .replace("ma_annum", "ma annum")
            .replace("unit = \"Gy\"", "unit = \"\"")
            .replace("unit = \"Gy / k_annum\"", "unit = \"\"");
        let mut inputs = io::inputs::parse_experiments(
            &fixture,
            1,
        )
        .expect("the checked-in input fixture should parse");
        let experiment = &mut inputs[0];
        experiment.cube.cube_height.value = 2.0;
        experiment.cube.cube_width.value = 2.0;
        experiment.cube.cube_length.value = 2.0;
        experiment.cube.hole_count = 1;
        experiment.cube.bandtail_count = 0;
        experiment.cube.periodic = false;
        experiment.traps[0].1.density.value = 1.0;
        experiment.traps[0].1.initial_fill.value = 0.25;
        experiment.transitions.ground_state_tunnelling = false;
        experiment.transitions.excited_state_tunnelling = false;
        experiment.transitions.vrh = false;
        experiment.transitions.ground_state_to_conduction_band = false;
        experiment.transitions.excited_state_to_conduction_band = false;
        experiment.transitions.filling = false;
        experiment.transitions.conduction_band_retrapping = false;
        experiment.transitions.ground_state_tunnelling_retrapping = false;
        experiment.transitions.excited_state_tunnelling_retrapping = false;
        experiment.transitions.filling_retrapping = false;
        inputs
    }

    #[test]
    fn new_builds_all_runtime_state_from_grouped_inputs() {
        let inputs = small_inputs();
        let simulation = MonteCarloSimulation::new(&inputs, 12, 1, 3)
            .expect("valid grouped inputs should build a simulation");

        assert_eq!(simulation.repetitions, 12);
        assert_eq!(simulation.experiments, 1);
        assert_eq!(simulation.cubes[0].trap_total, vec![8]);
        assert_eq!(simulation.cubes[0].hole_total, 8);
        assert_eq!(simulation.cubes[0].bandtail_total, 0);
        assert_eq!(simulation.time_temperature.len(), 1);
        assert_eq!(
            simulation.transitions[0].delocalised,
            DelocalisedRateEquation::None
        );
        assert_eq!(simulation.transitions[0].filling, FillingRateEquation::None);
    }

    #[test]
    fn generation_helpers_propagate_invalid_geometry_and_profiles() {
        let mut inputs = small_inputs();
        inputs[0].cube.cube_height.value = 0.0;
        let cube_error = MonteCarloSimulation::generate_cube_for_experiment(
            &inputs[0].cube,
            vec![1.0],
            1,
        )
        .unwrap_err();
        assert!(cube_error.contains("boundary x must be greater than zero"));

        let mut profile = inputs[0].time_temperature.clone();
        profile.temperatures.value.pop();
        assert!(MonteCarloSimulation::generate_time_temperature(&profile).is_err());
    }

    #[test]
    fn transition_generation_maps_flags_to_equations() {
        let mut inputs = small_inputs();
        inputs[0].transitions.ground_state_tunnelling = true;
        inputs[0].transitions.excited_state_tunnelling = false;
        inputs[0].transitions.ground_state_tunnelling_retrapping = false;
        inputs[0].transitions.excited_state_tunnelling_retrapping = true;
        inputs[0].transitions.ground_state_to_conduction_band = true;
        inputs[0].transitions.excited_state_to_conduction_band = false;
        inputs[0].transitions.filling = true;

        let transitions = MonteCarloSimulation::generate_transitions(&inputs[0].transitions).unwrap();
        assert_eq!(
            transitions.delocalised,
            DelocalisedRateEquation::Ground {
                re: DelocalisedRateEquationType::FirstOrder,
            },
        );
        assert_eq!(transitions.localised_recomb, LocalisedRateEquation::Ground);
        assert_eq!(transitions.localised_retrap, LocalisedRateEquation::Excited);
        assert_eq!(transitions.filling, FillingRateEquation::FirstOrder);
    }

    #[test]
    fn run_writes_every_experiment_and_repetition() {
        let inputs = small_inputs();
        let mut simulation = MonteCarloSimulation::new(&inputs, 3, 1, 1).unwrap();
        let output_directory = temporary_output_directory("simulation");
        fs::create_dir(&output_directory).unwrap();

        simulation.run_to_directory(&output_directory).unwrap();

        for repetition_index in 0..3 {
            let output_file = output_directory.join(format!(
                "experiment_results_0_{repetition_index}.bin.gz"
            ));
            assert!(output_file.is_file());
            fs::remove_file(output_file).unwrap();
        }
        fs::remove_dir(output_directory).unwrap();
    }

    #[test]
    fn final_ratios_are_grouped_by_experiment() {
        let inputs = small_inputs();
        let simulation = MonteCarloSimulation::new(&inputs, 3, 1, 1).unwrap();
        let ratios = simulation
            .run_to_final_ratio_only(simulation.time_temperature[0].clone(), 1)
            .unwrap();

        assert_eq!(ratios.len(), 1);
        assert_eq!(ratios[0].len(), 3);
        assert!(ratios[0].iter().all(|ratio| (ratio - 0.25).abs() < 1.0e-12));
    }
}
