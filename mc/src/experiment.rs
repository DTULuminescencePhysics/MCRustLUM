// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use crate::calculate_times::{RecordedEvent, run_standard};
use common::crystal::Cube;
use common::numeric::Float;
use common::rate_equation_selection::Transitions;
use common::time_temperature::TimeTemperature;
use common::trap_hole_band_tail::ElectronPlaces;
use io::inputs::SimulationInputs;
use io::outputs::create_monte_carlo_experiment_file;

use rand::Rng;
use serde::{Deserialize, Serialize};
use std::path::Path;
/// Holds the unique id of each trap. This currently is set to u16
/// limiting the number of traps to just over 65,000 which
/// should for now be sufficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct PlaceId(u16);

impl PlaceId {
    pub fn new(index: usize) -> Result<Self, String> {
        let index = u16::try_from(index).map_err(|_| "trap count exceeds u16::MAX".to_string())?;

        Ok(Self(index))
    }

    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Holds the PlaceIds for traps, holes or bandtail states.
/// Available places are kept at the front of the ids list and
/// currently unavailable places are at the back
/// [ available Ids | unavailable Ids ]
///                 ^
///           available_count
/// The two extremes of this are then
/// [ unavailable Ids ] and [ available Ids ]
/// ^                    |                  ^
/// available_count      |                  available_count
#[derive(Debug)]
pub struct PlaceAvailability {
    /// A permutation containing every PlaceIf exactly once.
    ids: Box<[PlaceId]>,
    /// positions[position_id] gives that place's current index in `ids`.
    positions: Box<[u16]>,
    /// ids[..available_count] are available.
    available_count: usize,
}

impl PlaceAvailability {
    /// Function that makes all places initially unavailable
    pub fn new(count: usize) -> Result<Self, String> {
        if count >= u16::MAX as usize {
            return Err("too many traps for u16 IDs".to_string());
        }

        let ids = (0..count)
            .map(PlaceId::new)
            .collect::<Result<Vec<_>, _>>()?
            .into_boxed_slice();

        let positions = (0..count)
            .map(|index| index as u16)
            .collect::<Vec<_>>()
            .into_boxed_slice();

        Ok(Self {
            ids,
            positions,
            available_count: 0,
        })
    }

    /// Randomly selects Ids to make them available at the beginning of an experiment
    pub fn set_initial_condition(count: usize, available_count: usize) -> Result<Self, String> {
        let mut places = PlaceAvailability::new(count)?;
        if available_count == 0 {
            return Ok(places);
        } else if available_count == count {
            places.mark_all_available();
            return Ok(places);
        } else {
            places.randomly_make_available(available_count)?;
            return Ok(places);
        }
    }
    pub fn randomly_make_available(&mut self, n: usize) -> Result<(), String> {
        let unavailable_count = self.ids.len() - self.available_count;
        if n > unavailable_count {
            return Err(format!(
                "cannot make {n} places available: only \
                {unavailable_count} places remain"
            ));
        }

        let first_new = self.available_count;
        let new_available_end = first_new + n;
        let mut rng = common::random::rng();
        for destination in first_new..new_available_end {
            let selected = rng.gen_range(destination..self.ids.len());
            self.swap_positions(destination, selected);
        }

        self.available_count = new_available_end;

        Ok(())
    }

    pub fn mark_all_occupied(&mut self) {
        self.available_count = 0;
    }

    pub fn mark_all_available(&mut self) {
        self.available_count = self.ids.len();
    }
    /// Gives Ids available for reaction
    /// i.e. an occupied trap or unoccupied hole
    pub fn available(&self) -> &[PlaceId] {
        &self.ids[..self.available_count]
    }
    /// Gives Ids not currently available for reaction
    /// i.e. an unoccupied trap or occupied hole
    pub fn unavailable(&self) -> &[PlaceId] {
        &self.ids[self.available_count..]
    }

    /// Checks if a given PlaceId is available to the program
    pub fn is_available(&self, place: PlaceId) -> bool {
        if (self.positions[place.index()] as usize) < self.available_count {
            return true;
        } else {
            return false;
        }
    }
    /// Returns the availability count
    pub fn available_count(&self) -> usize {
        self.available_count
    }

    pub fn fill_ratio(&self) -> Float {
        self.available_count() as Float / self.ids.len() as Float
    }

    /// Swaps the two entries
    fn swap_positions(&mut self, first: usize, second: usize) {
        if first == second {
            return;
        }

        self.ids.swap(first, second);

        let first_id = self.ids[first];
        let second_id = self.ids[second];

        self.positions[first_id.index()] = first as u16;
        self.positions[second_id.index()] = second as u16;
    }

    /// To make a PlaceId available it needs to be swapped with the first
    /// unavailable id and the available count increased
    /// [ available Ids | A, C, ... unavailable Ids... B, ... ]
    ///                   ^                            ^
    ///            first unavailable Id             Id to move
    /// [ available Ids | B, C, ... unavailable Ids... A, ... ]
    ///                   ^                            ^
    ///               moved Id                former first unavailable Id
    /// [ available Ids B | C, ... unavailable Ids... A, ... ]
    ///                 ^   ^
    ///          moved Id   new first unavailable Id
    pub fn make_available(&mut self, place: PlaceId) -> bool {
        let current = self.positions[place.index()] as usize;

        if current < self.available_count {
            return false; // Already available
        }

        let first_occupied = self.available_count;
        self.swap_positions(current, first_occupied);
        self.available_count += 1;

        true
    }
    /// To make a PlaceId unavailable it needs to be swapped with the last
    /// available id and the available count decreased
    /// [ B, ... available Ids ... C, A | unavailable Ids ]
    ///   ^                           ^
    ///   Id to move         last available Id
    /// [ A, ... available Ids ... C, B | unavailable Ids ]
    ///   ^                           ^
    /// former last available Id   moved Id    
    /// [ A, ... available Ids ... C | B, ... unavailable Ids ]
    ///                            ^   ^
    ///        new last available Id   moved Id              
    pub fn make_unavailable(&mut self, trap: PlaceId) -> bool {
        let current = self.positions[trap.index()] as usize;

        if current >= self.available_count {
            return false; // Already occupied
        }

        let last_available = self.available_count - 1;
        self.swap_positions(current, last_available);
        self.available_count = last_available;

        true
    }
}

/// The trap parameters
#[derive(Debug, Clone, Copy)]
pub struct TrapParameters {
    pub(crate) excited_energy_gap: Float,
    pub(crate) s_frequency_e: Float,
    pub(crate) s_frequency_g: Float,
    pub(crate) e_cb_ground: Float,
    pub(crate) e_cb_excited: Float,
    pub(crate) de_frequency_ground: Float,
    pub(crate) de_frequency_excited: Float,
    pub(crate) lo_frequency_ground: Float,
    pub(crate) lo_frequency_excited: Float,
    pub(crate) alpha_ground: Float,
    pub(crate) alpha_excited: Float,
    pub(crate) delocalised_mu: Float,
    pub(crate) retrap_ratio: Float,
}
impl TrapParameters {
    pub fn new(
        excited_energy_gap: Float,
        s_frequency_e: Float,
        s_frequency_g: Float,
        e_cb_ground: Float,
        e_cb_excited: Float,
        de_frequency_ground: Float,
        de_frequency_excited: Float,
        lo_frequency_ground: Float,
        lo_frequency_excited: Float,
        alpha_ground: Float,
        alpha_excited: Float,
        delocalised_mu: Float,
        retrap_ratio: Float,
    ) -> Self {
        Self {
            excited_energy_gap,
            s_frequency_e,
            s_frequency_g,
            e_cb_ground,
            e_cb_excited,
            de_frequency_ground,
            de_frequency_excited,
            lo_frequency_ground,
            lo_frequency_excited,
            alpha_ground,
            alpha_excited,
            delocalised_mu,
            retrap_ratio,
        }
    }
}

pub enum TrapParameterLayout {
    /// Every trap uses exactly the same parameters.
    Uniform(TrapParameters),

    /// Every trap has its own parameters.
    ///
    /// parameters[trap_id]
    Direct(Box<[TrapParameters]>),

    /// Traps reference a table of shared parameter records.
    ///
    /// Useful for families and mixtures of shared/individual parameters.
    Indexed {
        records: Box<[TrapParameters]>,
        by_trap: Box<[PlaceId]>,
    },
}

impl TrapParameterLayout {
    pub fn new_uniform(inputs: &SimulationInputs) -> Self {
        let parameters = TrapParameters::new(
            inputs.trap_energies.e_loc[0],
            inputs.trap_energies.s_frequency_e[0],
            inputs.trap_energies.s_frequency_g[0],
            inputs.trap_energies.e_cb[0],
            inputs.trap_energies.e_cb[0] - inputs.trap_energies.e_loc[0],
            inputs.delocalised.s_gs[0],
            inputs.delocalised.s_es[0],
            inputs.localised.b_gs[0],
            inputs.localised.b_es[0],
            inputs.localised.alpha_gs[0],
            inputs.localised.alpha_es[0],
            inputs.delocalised.mu[0],
            inputs.delocalised.retrap_ratio[0],
        );
        TrapParameterLayout::uniform(parameters)
    }

    pub fn uniform(parameters: TrapParameters) -> Self {
        Self::Uniform(parameters)
    }

    pub fn direct(parameters: Vec<TrapParameters>, trap_count: usize) -> Result<Self, String> {
        if parameters.len() != trap_count {
            return Err(format!(
                "expected {trap_count} trap parameter records, found {}",
                parameters.len(),
            ));
        }

        Ok(Self::Direct(parameters.into_boxed_slice()))
    }

    pub fn indexed(
        records: Vec<TrapParameters>,
        assignments: Vec<PlaceId>,
        trap_count: usize,
    ) -> Result<Self, String> {
        if records.is_empty() {
            return Err("at least one parameter record is required".to_string());
        }

        if assignments.len() != trap_count {
            return Err(format!(
                "expected {trap_count} parameter assignments, found {}",
                assignments.len(),
            ));
        }

        for (trap_index, parameter_id) in assignments.iter().enumerate() {
            if parameter_id.index() >= records.len() {
                return Err(format!(
                    "trap {trap_index} references missing parameter {}",
                    parameter_id.index(),
                ));
            }
        }

        Ok(Self::Indexed {
            records: records.into_boxed_slice(),
            by_trap: assignments.into_boxed_slice(),
        })
    }

    pub fn get(&self, trap: PlaceId) -> &TrapParameters {
        match self {
            Self::Uniform(parameters) => parameters,

            Self::Direct(parameters) => &parameters[trap.index()],

            Self::Indexed { records, by_trap } => &records[by_trap[trap.index()].index()],
        }
    }
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
    },
    /// Contains parts for experiment with bandtails
    WithBandtail {
        places: ElectronPlaces,
        trap_places: PlaceAvailability,
        hole_places: PlaceAvailability,
        bandtail_places: PlaceAvailability,
        trap_parameters: TrapParameterLayout,
        time_temperature: TimeTemperature,
    },
}

impl MCExperiment {
    pub fn initialise(
        cube: &Cube,
        inputs: &SimulationInputs,
        time_temperature: TimeTemperature,
    ) -> Result<Self, String> {
        let places = ElectronPlaces::random_from_cube(cube)?;
        let trap_places = PlaceAvailability::new(cube.trap_total)?;
        let hole_places = PlaceAvailability::new(cube.hole_total)?;
        let trap_parameters = TrapParameterLayout::new_uniform(inputs);

        if cube.bandtail_total == 0 {
            Ok(Self::Standard {
                places,
                trap_places,
                hole_places,
                trap_parameters,
                time_temperature,
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
    use crate::calculate_times::Event;
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
        let mut experiment = MCExperiment::initialise(&cube, &inputs, profile).unwrap();
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
        assert_eq!(batches[0].len(), 1);
        assert_eq!(batches[0][0].event, Event::None);
        assert_eq!(batches[0][0].time, 1.0);
    }
}
