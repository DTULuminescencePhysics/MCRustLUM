// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use crate::calculate_times::run_standard;
use common::charge_transfer::RecordedEvent;
use common::crystal::Cube;
use common::numeric::Float;
use common::place_ids::PlaceAvailability;
use common::rate_equation_selection::Transitions;
use common::time_temperature::TimeTemperature;
use common::trap_hole_band_tail::{ElectronPlaces, TrapParameterLayout};
use io::inputs::SimulationInputs;
use io::outputs::create_monte_carlo_experiment_file;
use std::path::Path;
use rand::rngs::StdRng;
use common::random::get_std_rng_for_rep;

pub fn experiment_value(
    values: &[Float],
    experiment_index: usize,
    field: &str,
) -> Result<Float, String> {
    match values {
        [] => Err(format!("{field} must contain at least one value")),
        [value] => Ok(*value),
        _ => values.get(experiment_index).copied().ok_or_else(|| {
            format!(
                "{field} contains {} values, but experiment {} was requested",
                values.len(),
                experiment_index + 1,
            )
        }),
    }
}

pub fn new_uniform_trap_layout(
    inputs: &SimulationInputs,
    exp: &usize,
) -> Result<TrapParameterLayout, String> {
    Ok(TrapParameterLayout::new_uniform(
        experiment_value(&inputs.trap_energies.e_loc, *exp, "trap_energies.e_loc")?,
        experiment_value(
            &inputs.trap_energies.s_frequency_e,
            *exp,
            "trap_energies.s_frequency_e",
        )?,
        experiment_value(
            &inputs.trap_energies.s_frequency_g,
            *exp,
            "trap_energies.s_frequency_g",
        )?,
        experiment_value(&inputs.trap_energies.e_cb, *exp, "trap_energies.e_cb")?,
        experiment_value(&inputs.delocalised.s_gs, *exp, "delocalised.s_gs")?,
        experiment_value(&inputs.delocalised.s_es, *exp, "delocalised.s_es")?,
        experiment_value(&inputs.localised.b_gs, *exp, "localised.b_gs")?,
        experiment_value(&inputs.localised.b_es, *exp, "localised.b_es")?,
        experiment_value(&inputs.localised.alpha_gs, *exp, "localised.alpha_gs")?,
        experiment_value(&inputs.localised.alpha_es, *exp, "localised.alpha_es")?,
        experiment_value(&inputs.delocalised.mu, *exp, "delocalised.mu")?,
        experiment_value(
            &inputs.delocalised.retrap_ratio,
            *exp,
            "delocalised.retrap_ratio",
        )?,
    ))
}

/// Holds
pub enum MCExperiment {
    /// Contains experiment parts for standard run
    Standard {
        places: ElectronPlaces,
        trap_places: PlaceAvailability,
        hole_places: PlaceAvailability,
        trap_parameters: TrapParameterLayout,
        time_temperature: TimeTemperature,
        rng: StdRng,
    },
    /// Contains parts for experiment with bandtails
    WithBandtail {
        places: ElectronPlaces,
        trap_places: PlaceAvailability,
        hole_places: PlaceAvailability,
        bandtail_places: PlaceAvailability,
        trap_parameters: TrapParameterLayout,
        time_temperature: TimeTemperature,
        rng: StdRng,
    },
}

impl MCExperiment {
    pub fn initialise(
        cube: &Cube,
        inputs: &SimulationInputs,
        trap_available: &usize,
        hole_available: &usize,
        time_temperature: TimeTemperature,
        exp: &usize,
        rep: &usize,
        
    ) -> Result<Self, String> {

        let mut rng = get_std_rng_for_rep(*rep);
        let places = ElectronPlaces::random_from_cube(cube, &mut rng)?;
        let trap_places =
            PlaceAvailability::set_initial_condition(cube.trap_total, *trap_available, &mut rng)?;
        let hole_places =
            PlaceAvailability::set_initial_condition(cube.hole_total, *hole_available, &mut rng)?;
        let trap_parameters = new_uniform_trap_layout(inputs, exp)?;

        if cube.bandtail_total == 0 {
            Ok(Self::Standard {
                places,
                trap_places,
                hole_places,
                trap_parameters,
                time_temperature,
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
                rng
            })
        }
    }

    pub fn run(
        &mut self,
        cube: &Cube,
        inputs: &SimulationInputs,
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
                    time_temperature,
                    cube,
                    inputs,
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::charge_transfer::Event;
    use common::constants::temperature::TemperatureUnit;
    use common::constants::time::TimeUnit;
    use io::outputs::read_all_batches;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_output_path(label: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_{label}_{}_{}.bin.gz",
            std::process::id(),
            unique,
        ))
    }

    #[test]
    fn standard_run_dispatches_to_run_standard() {
        let cube = Cube::new(1.0, 1.0, 1.0, 1, 1, 0, false).unwrap();
        let inputs = SimulationInputs::default();
        let profile = TimeTemperature::new(
            vec![0.0, 1.0],
            vec![20.0, 20.0],
            TimeUnit::Second,
            TemperatureUnit::Celsius,
        )
        .unwrap();
        let transitions = Transitions::from_bool(
            false, false, false, false, false, "first", false, false, false, false,
        )
        .unwrap();
        let mut experiment = MCExperiment::initialise(&cube, &inputs, &0, &0, profile, &0, &0).unwrap();
        let output_file = temporary_output_path("standard_experiment");

        experiment
            .run(&cube, &inputs, &transitions, &output_file, &2)
            .unwrap();
        let batches = read_all_batches::<RecordedEvent>(&output_file)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        fs::remove_file(output_file).unwrap();

        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].len(), 2);
        assert_eq!(batches[0][0].event, Event::None);
        assert_eq!(batches[0][0].time, 0.0);
        assert_eq!(batches[0][1].event, Event::None);
        assert_eq!(batches[0][1].time, 1.0);
    }

    #[test]
    fn singleton_parameters_are_shared_between_experiments() {
        let inputs = SimulationInputs::default();

        assert!(new_uniform_trap_layout(&inputs, &1).is_ok());
    }

    #[test]
    fn short_parameter_vectors_return_an_error_instead_of_panicking() {
        let mut inputs = SimulationInputs::default();
        inputs.trap_energies.e_loc = vec![1.0, 2.0];

        let error = new_uniform_trap_layout(&inputs, &2)
            .err()
            .expect("a short parameter vector should be rejected");

        assert_eq!(
            error,
            "trap_energies.e_loc contains 2 values, but experiment 3 was requested"
        );
    }
}
