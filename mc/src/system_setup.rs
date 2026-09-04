// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Creation and per-iteration reset of Monte Carlo state.

use crate::experiment::{MCExperiment, experiment_value};
use common::crystal::Cube;
use common::numeric::Float;
use common::rate_equation_selection::Transitions;
use common::time_temperature::TimeTemperature;
use io::inputs::{CubeSpecification, SimulationInputs, TimeTempSpecification};
use std::path::Path;

/// State shared by the repetitions and experiments in one Monte Carlo run.
///
/// Constructing this type validates the supplied geometry and temperature
/// profile and randomly places the crystal sites.
#[derive(Debug, Clone)]
pub struct MonteCarloSimulation {
    /// Original grouped configuration, retained for regeneration and reset.
    pub inputs: SimulationInputs,
    /// Current spatial realisation of traps, holes, and bandtail states.
    pub cube: Cube,
    /// Current time and temperature state.
    pub time_temperature: TimeTemperature,
    /// Enabled physical pathways and their rate equations.
    pub transitions: Transitions,
    /// Number of repetitions performed within each experiment.
    pub repetions: usize,
    /// Number of independent experiments to run.
    pub experiments: usize,
}

impl MonteCarloSimulation {
    /// Build all runtime state from a grouped input configuration.
    ///
    /// ```no_run
    /// # fn main() -> Result<(), String> {
    /// let simulation = mc::system_setup::MonteCarloSimulation::new(
    ///     io::default_inputs(),
    ///     10,
    ///     1,
    /// )?;
    /// assert_eq!(simulation.repetions, 10);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        inputs: SimulationInputs,
        repetions: usize,
        experiments: usize,
    ) -> Result<Self, String> {
        let cube = MonteCarloSimulation::generate_cube(&inputs.cube)?;
        let time_temperature =
            MonteCarloSimulation::generate_time_temperature(&inputs.time_temperature)?;
        let transitions = MonteCarloSimulation::generate_transitions(&inputs)?;

        Ok(Self {
            inputs,
            cube,
            time_temperature,
            transitions,
            repetions,
            experiments,
        })
    }

    /// Generate a cube with new random site coordinates.
    ///
    /// The trap count is derived from the configured density and volume. Hole
    /// and bandtail counts are interpreted as ratios per trap.
    pub fn generate_cube(inputs: &CubeSpecification) -> Result<Cube, String> {
        Cube::new_from_density(
            inputs.x,
            inputs.y,
            inputs.z,
            inputs.density,
            inputs.hole_count,
            inputs.bandtail_count,
            inputs.periodic,
        )
    }
    /// Validate and generate the piecewise-linear time/temperature profile.
    pub fn generate_time_temperature(
        inputs: &TimeTempSpecification,
    ) -> Result<TimeTemperature, String> {
        TimeTemperature::new(
            inputs.times.clone(),
            inputs.temperatures.clone(),
            inputs.time_unit,
            inputs.temp_unit,
        )
    }
    /// Generate the transition selection used by the simulation.
    ///
    /// Delocalised release currently uses first-order kinetics. Boolean flags
    /// from the localised, delocalised, and filling input groups determine
    /// which pathways are active.
    pub fn generate_transitions(inputs: &SimulationInputs) -> Result<Transitions, String> {
        Transitions::from_bool(
            inputs.localised.gs_tun,
            inputs.localised.es_tun,
            inputs.delocalised.gs_cb,
            inputs.delocalised.es_cb,
            inputs.filling.fill,
            "first",
            inputs.localised.gs_retrap,
            inputs.localised.es_retrap,
            inputs.delocalised.retrap,
            inputs.filling.cmbn_whn_fll,
        )
    }

    /// Restore the time/temperature profile to its first control point.
    pub fn reset_time_temperature(&mut self) {
        self.time_temperature.reset();
    }

    /// Run every configured experiment and repetition.
    ///
    /// Every repetition receives a newly generated random spatial realization
    /// together with a reset copy of the time/temperature profile.
    pub fn run(&self) -> Result<(), String> {
        self.run_to_directory(Path::new("tmp"))
    }

    /// Run every configured experiment and repetition, placing temporary
    /// result files in `output_directory`.
    pub fn run_to_directory(&self, output_directory: impl AsRef<Path>) -> Result<(), String> {
        let output_directory = output_directory.as_ref();
        let batch_capacity: usize = 100;
        for experiment_index in 0..self.experiments {
            let trap_available = (experiment_value(
                &self.inputs.initial_conditions.trap_available,
                experiment_index,
                "initial_conditions.trap_available",
            )? * self.cube.trap_total as Float) as usize;
            let hole_available = (experiment_value(
                &self.inputs.initial_conditions.hole_available,
                experiment_index,
                "initial_conditions.hole_available",
            )? * self.cube.hole_total as Float) as usize;
            let experiment_offset = experiment_index*self.repetions; 
            for repetition_index in 0..self.repetions {
                let temp_file_path = output_directory.join(format!(
                    "experiment_results_{}_{}.bin.gz",
                    experiment_index, repetition_index
                ));
                let mut time_temperature = self.time_temperature.clone();
                time_temperature.reset();
                let randrep = &experiment_offset + repetition_index;
                let mut experiment = MCExperiment::initialise(
                    &self.cube,
                    &self.inputs,
                    &trap_available,
                    &hole_available,
                    time_temperature,
                    &experiment_index,
                    &randrep,

                )
                .map_err(|error| {
                    format!(
                        "could not initialise experiment {} of {}, repetition {} of {}: {error}",
                        experiment_index + 1,
                        self.experiments,
                        repetition_index + 1,
                        self.repetions,
                    )
                })?;
                experiment
                    .run(
                        &self.cube,
                        &self.inputs,
                        &self.transitions,
                        &temp_file_path,
                        &batch_capacity,
                    )
                    .map_err(|error| {
                        format!(
                            "experiment {} of {}, repetition {} of {} failed: {error}",
                            experiment_index + 1,
                            self.experiments,
                            repetition_index + 1,
                            self.repetions,
                        )
                    })?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::charge_transfer::{Event, RecordedEvent};
    use common::rate_equation_selection::{
        DelocalisedRateEquation, DelocalisedRateEquationType, FillingRateEquation,
        LocalisedRateEquation,
    };
    use io::outputs::read_all_batches;
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

    fn small_inputs() -> SimulationInputs {
        let mut inputs = io::default_inputs();
        inputs.cube.x = 2.0;
        inputs.cube.y = 1.0;
        inputs.cube.z = 1.0;
        inputs.cube.density = 1.5;
        inputs.cube.hole_count = 2;
        inputs.cube.bandtail_count = 1;
        inputs.cube.periodic = false;
        inputs.time_temperature.times = vec![0.0, 10.0];
        inputs.time_temperature.temperatures = vec![20.0, 30.0];
        inputs
    }

    #[test]
    fn new_builds_all_runtime_state_from_grouped_inputs() {
        let simulation = MonteCarloSimulation::new(small_inputs(), 12, 3)
            .expect("valid grouped inputs should build a simulation");

        assert_eq!(simulation.repetions, 12);
        assert_eq!(simulation.experiments, 3);
        assert_eq!(simulation.cube.trap_total, 3);
        assert_eq!(simulation.cube.hole_total, 6);
        assert_eq!(simulation.cube.bandtail_total, 3);
        assert!(!simulation.inputs.cube.periodic);
        assert!((simulation.time_temperature.current_temperature() - 293.15).abs() < 1.0e-12);
        assert!(matches!(
            simulation.transitions,
            Transitions::CbRetrapping { .. }
        ));
    }

    #[test]
    fn generation_helpers_propagate_invalid_geometry_and_profiles() {
        let mut cube = CubeSpecification::default();
        cube.x = 0.0;
        let cube_error = MonteCarloSimulation::generate_cube(&cube).unwrap_err();
        assert!(cube_error.contains("boundary x must be greater than zero"));

        let mut profile = TimeTempSpecification::default();
        profile.temperatures = vec![20.0];
        let profile_error = MonteCarloSimulation::generate_time_temperature(&profile).unwrap_err();
        assert_eq!(
            profile_error,
            "times and temperatures must have the same length"
        );
    }

    #[test]
    fn transition_generation_maps_flags_to_equations_and_outer_variant() {
        let mut inputs = small_inputs();
        inputs.localised.gs_tun = true;
        inputs.localised.es_tun = false;
        inputs.localised.gs_retrap = false;
        inputs.localised.es_retrap = true;
        inputs.delocalised.gs_cb = true;
        inputs.delocalised.es_cb = false;
        inputs.delocalised.retrap = true;
        inputs.filling.fill = true;
        inputs.filling.cmbn_whn_fll = true;

        let transitions = MonteCarloSimulation::generate_transitions(&inputs).unwrap();
        match transitions {
            Transitions::FillCbRetrapping { transitions } => {
                assert_eq!(
                    transitions.delocalised,
                    DelocalisedRateEquation::Ground {
                        re: DelocalisedRateEquationType::FirstOrder,
                    },
                );
                assert_eq!(transitions.localised_recomb, LocalisedRateEquation::Ground);
                assert_eq!(transitions.localised_retrap, LocalisedRateEquation::Excited);
                assert_eq!(transitions.filling, FillingRateEquation::Basic);
            }
            other => panic!("expected both retrapping flags, got {other:?}"),
        }
    }

    #[test]
    fn run_writes_every_experiment_and_repetition() {
        let mut inputs = small_inputs();
        inputs.cube.x = 1.0;
        inputs.cube.y = 1.0;
        inputs.cube.z = 1.0;
        inputs.cube.density = 1.0;
        inputs.cube.hole_count = 1;
        inputs.cube.bandtail_count = 0;
        inputs.time_temperature.times = vec![0.0, 1.0];
        inputs.time_temperature.temperatures = vec![20.0, 20.0];
        inputs.initial_conditions.trap_available = vec![1.0, 1.0];
        inputs.initial_conditions.hole_available = vec![1.0, 1.0];
        inputs.localised.gs_tun = false;
        inputs.localised.es_tun = false;
        inputs.localised.gs_retrap = false;
        inputs.localised.es_retrap = false;
        inputs.delocalised.gs_cb = false;
        inputs.delocalised.es_cb = false;
        inputs.delocalised.retrap = false;
        inputs.filling.fill = false;
        inputs.filling.cmbn_whn_fll = false;

        let simulation = MonteCarloSimulation::new(inputs, 3, 2).unwrap();
        let output_directory = temporary_output_directory("simulation");
        fs::create_dir(&output_directory).unwrap();

        simulation.run_to_directory(&output_directory).unwrap();

        for experiment_index in 0..2 {
            for repetition_index in 0..3 {
                let output_file = output_directory.join(format!(
                    "experiment_results_{experiment_index}_{repetition_index}.bin.gz"
                ));
                let batches = read_all_batches::<RecordedEvent>(&output_file)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();

                assert_eq!(batches.len(), 1);
                assert_eq!(batches[0].len(), 2);
                assert_eq!(batches[0][0].event, Event::None);
                assert_eq!(batches[0][1].event, Event::None);
                fs::remove_file(output_file).unwrap();
            }
        }
        fs::remove_dir(output_directory).unwrap();
    }
}
