//! Creation and per-iteration reset of Monte Carlo state.

use common::crystal:: Cube;
use common::numeric::Float;
use common::time_temperature::TimeTemperature;
use io::inputs::{SimulationInputs, CubeSpecification, TimeTempSpecification};

/// The immutable information needed to generate a fresh random [`Cube`].
///
/// Keeping this separately from the active cube means an iteration may add or
/// remove sites without changing the number of sites generated for the next
/// iteration.


#[derive(Debug, Clone)]
pub struct MonteCarloSimulation{
    pub inputs: SimulationInputs,
    pub cube: Cube, 
    pub time_temperature: TimeTemperature,
    pub experiments: usize,
    pub repetions: usize,

}

impl MonteCarloSimulation {
    /// Generate a cube with new random site coordinates.
    pub fn generate_cube(inputs: CubeSpecification) -> Cube {
        Cube::new_random_from_density(
            inputs.x,
            inputs.y,
            inputs.z,
            inputs.density,
            inputs.hole_count,
            inputs.bandtail_count,
            inputs.periodic,
        )
    }

    pub fn generate_time_temperature(inputs: TimeTempSpecification) -> Result<TimeTemperature, String>  {
        TimeTemperature::new(
            inputs.times,
            inputs.temperatures,
            inputs.time_unit,
            inputs.temp_unit,
        )
    }

    /// Replace the current cube with a new random spatial realization.
    pub fn regenerate_cube(&mut self) {
        self.cube.randomise_positions();
    }

    /// Restore the time/temperature profile to its first control point.
    pub fn reset_time_temperature(&mut self) {
        self.time_temperature.reset();
    }

}




