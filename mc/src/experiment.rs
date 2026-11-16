// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Construction and execution of one Monte Carlo repetition.
//!
//! An experiment owns its independently randomized site coordinates,
//! occupancy partitions, temperature-profile state, and random-number stream.
//! Trap-family parameters may vary between experiment indices while remaining
//! uniform across all traps within one current experiment.

use crate::calculate_times::{run_standard, run_final_ratio_only};
use common::charge_transfer::RecordedEvent;
use common::crystal::Cube;
use common::numeric::Float;
use common::place_ids::PlaceAvailability;
use common::rate_equation_selection::Transitions;
use common::time_temperature::TimeTemperature;
use common::trap_hole_band_tail::{TrapGenerator, ElectronPlaces, TrapParameterLayout, ReTrapParameterLayout};
use common::rate_equation_inputs::FillingTransitionInputs;
use io::outputs::create_monte_carlo_experiment_file;
use std::path::Path;
use rand::rngs::StdRng;
use common::random::get_std_rng_for_rep;

/// Mutable state for one independently randomized Monte Carlo repetition.
///
/// The variant records whether the spatial model includes shallow band-tail
/// sites. Only [`MCExperiment::Standard`] currently has a kinetic runner.
pub enum MCExperiment {
    /// A trap-and-hole model without band-tail transport.
    Standard {
        /// Randomized trap and hole coordinates.
        places: ElectronPlaces,
        /// Partition whose available traps are electron occupied.
        trap_places: PlaceAvailability,
        /// Partition whose available sites contain recombining holes.
        hole_places: PlaceAvailability,
        /// Physical parameter assignment for every trap.
        trap_parameters: TrapParameterLayout,
        /// Independent mutable copy of the experimental profile.
        time_temperature: TimeTemperature,
        /// Retrapping parameters
        retrapping_parameters: ReTrapParameterLayout,
        /// Filling parameters
        filling_parameters: Option<FillingTransitionInputs>,
        /// Repetition-specific random-number generator.
        rng: StdRng,
    },
    /// A spatial model that also contains band-tail states.
    WithBandtail {
        /// Randomized trap, hole, and band-tail coordinates.
        places: ElectronPlaces,
        /// Partition whose available traps are electron occupied.
        trap_places: PlaceAvailability,
        /// Partition whose available sites contain recombining holes.
        hole_places: PlaceAvailability,
        /// Occupancy partition for band-tail states.
        bandtail_places: PlaceAvailability,
        /// Physical parameter assignment for every trap.
        trap_parameters: TrapParameterLayout,
        /// Independent mutable copy of the experimental profile.
        time_temperature: TimeTemperature,
        /// Retrapping parameters
        retrapping_parameters: ReTrapParameterLayout,
        /// Filling parameters
        filling_parameters: Option<FillingTransitionInputs>,
        /// Repetition-specific random-number generator.
        rng: StdRng,
    },
}

impl MCExperiment {
    /// Randomize site positions and initial populations for one repetition.
    ///
    /// `trap_available` is the requested number of initially occupied traps;
    /// `hole_available` is the number of active hole centres. The `rep` index
    /// derives a deterministic random seed, while `exp` selects any
    /// experiment-dependent physical parameter values.
    pub fn initialise(
        cube: &Cube,
        trap_available: &Vec<usize>,
        hole_available: &usize,
        trap_generator: &Vec<TrapGenerator>,
        retrapping_parameters: ReTrapParameterLayout,
        filling: &mut Option<FillingTransitionInputs>,
        time_temperature: TimeTemperature,
        rep: &usize,
        
    ) -> Result<Self, String> {

        let mut rng = get_std_rng_for_rep(*rep);

        let trap_parameters = TrapParameterLayout::new(
            trap_generator, 
            &cube.trap_total,
            &mut rng )?;
        
        let places = ElectronPlaces::random_from_cube(cube, &mut rng)?;     

        let trap_places = PlaceAvailability::set_initial_condition(
            &cube.trap_total,
            trap_available, 
            &mut rng)?;
        let hole_places = PlaceAvailability::set_initial_condition(
            &vec![cube.hole_total], 
            &vec![*hole_available], 
            &mut rng)?;
        
        let filling_parameters = match filling {
            Some(filling) => {
                filling.occupied_population = trap_available.iter().sum::<usize>() as Float;
                Some(*filling)
            }
            None => None
        };
        
       
        if cube.bandtail_total == 0 {
            Ok(Self::Standard {
                places,
                trap_places,
                hole_places,
                trap_parameters,
                time_temperature,
                retrapping_parameters,
                filling_parameters,
                rng
            })
        } else {
            let bandtail_places = PlaceAvailability::new(cube.bandtail_total)?;
            Ok(Self::WithBandtail {
                places,
                trap_places,
                hole_places,
                bandtail_places,
                trap_parameters,
                time_temperature,
                retrapping_parameters,
                filling_parameters,
                rng
            })
        }
    }

    /// Execute this repetition and stream its records to `output_file`.
    ///
    /// The output is truncated before a standard run and then written in
    /// bounded batches. Band-tail kinetics currently return an explicit error
    /// without creating output.
    pub fn run(
        &mut self,
        cube: &Cube,
        transitions: &Transitions,
        output_file: &Path,
        batch_capacity: &usize,
    ) -> Result<(), String> {
        match self {
            Self::Standard {
                places,
                trap_places,
                hole_places,
                trap_parameters,
                time_temperature,
                filling_parameters,
                retrapping_parameters,
                rng
            } => {
                create_monte_carlo_experiment_file(output_file)
                    .map_err(|error| error.to_string())?;
                let results: Vec<RecordedEvent> = Vec::with_capacity(*batch_capacity);
                run_standard(
                    places,
                    trap_places,
                    hole_places,
                    trap_parameters,
                    retrapping_parameters,
                    filling_parameters,
                    time_temperature,
                    cube,
                    transitions,
                    output_file,
                    results,
                    rng,
                )
            }
            Self::WithBandtail { .. } => {
                Err("bandtail kinetic Monte Carlo runs are not implemented yet".to_string())
            }
        }
    }

    /// Execute this repetition but only the final ratio is stored.
    ///
    /// Band-tail kinetics currently return an explicit error
    /// without creating output.
    pub fn run_no_output_file(
        &mut self,
        cube: &Cube,
        transitions: &Transitions,
    ) -> Result<Float, String> {
        match self {
            Self::Standard {
                places,
                trap_places,
                hole_places,
                trap_parameters,
                time_temperature,
                filling_parameters,
                retrapping_parameters,
                rng
            } => {
                run_final_ratio_only(
                    places,
                    trap_places,
                    hole_places,
                    trap_parameters,
                    retrapping_parameters,
                    filling_parameters,
                    time_temperature,
                    cube,
                    transitions,
                    rng,
                )
            }
            Self::WithBandtail { .. } => {
                Err("bandtail kinetic Monte Carlo runs are not implemented yet".to_string())
            }
        }
    }



}

#[cfg(test)]
mod tests {
    use super::*;
    use common::constants::time::TimeUnit;
    use common::crystal::Cube;
    use common::rate_equation_inputs::FillingTransitionInputs;
    use common::trap_hole_band_tail::{ReTrapParameterLayout, TrapGenerator};

    #[test]
    /// Verify that initialise uses current runtime inputs and updates filling population.
    fn initialise_uses_current_runtime_inputs_and_updates_filling_population() {
        let cube = Cube::new_from_density(1.0, 1.0, 1.0, vec![1.0], 1, 0, false, 1).unwrap();
        let trap_generators = vec![TrapGenerator::new_no_distribution(
            "test".into(), 0.1, 1.0, 1.0, 1.0, 0.9,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        )];
        let profile = TimeTemperature::new(
            vec![0.0, 1.0],
            vec![293.15, 293.15],
            TimeUnit::Second,
        )
        .unwrap();
        let mut filling = Some(FillingTransitionInputs::get_inputs(10.0, 1.0, 0, 1).unwrap());
        let trap_available = vec![1];

        let experiment = MCExperiment::initialise(
            &cube,
            &trap_available,
            &1,
            &trap_generators,
            ReTrapParameterLayout::None,
            &mut filling,
            profile,
            &0,
        )
        .unwrap();

        assert!(matches!(experiment, MCExperiment::Standard { .. }));
        assert_eq!(filling.unwrap().occupied_population, 1.0);
    }
}
