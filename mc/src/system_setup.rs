//! Creation and per-iteration reset of Monte Carlo state.

use common::crystal::Cube;
use common::time_temperature::TimeTemperature;
use common::rate_equation_selection::Transitions;
use io::inputs::{SimulationInputs, CubeSpecification, TimeTempSpecification};
/// State shared by the repetitions and experiments in one Monte Carlo run.
///
/// Constructing this type validates the supplied geometry and temperature
/// profile and randomly places the crystal sites.
#[derive(Debug, Clone)]
pub struct MonteCarloSimulation{
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
       Transitions::from_bool(inputs.localised.gs_tun, inputs.localised.es_tun,
                              inputs.delocalised.gs_cb, inputs.delocalised.es_cb,
                              inputs.filling.fill, "first",
                              inputs.localised.gs_retrap, inputs.localised.es_retrap,
                              inputs.delocalised.retrap,
                              inputs.filling.cmbn_whn_fll)

    }

    /// Replace the current cube with a new random spatial realisation.
    ///
    /// Counts, dimensions, and boundary conditions remain unchanged; cached
    /// distance matrices are recalculated for the new positions.
    pub fn regenerate_cube(&mut self) -> Result<(), String> {
        self.cube.randomise_positions()
    }

    /// Restore the time/temperature profile to its first control point.
    pub fn reset_time_temperature(&mut self) {
        self.time_temperature.reset();
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use common::rate_equation_selection::{
        DelocalisedRateEquation, DelocalisedRateEquationType, FillingRateEquation,
        LocalisedRateEquation,
    };

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
        assert_eq!(simulation.cube.places.traps.len(), 3);
        assert_eq!(simulation.cube.places.holes.len(), 6);
        assert_eq!(simulation.cube.places.bandtails.len(), 3);
        assert!(!simulation.inputs.cube.periodic);
        assert!((simulation.time_temperature.current_temperature() - 293.15).abs() < 1.0e-12);
        assert!(matches!(simulation.transitions, Transitions::CbRetrapping { .. }));
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
        assert_eq!(profile_error, "times and temperatures must have the same length");
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
    fn regeneration_and_reset_preserve_configuration_but_replace_state() {
        let mut simulation = MonteCarloSimulation::new(small_inputs(), 1, 1).unwrap();
        let original: Vec<_> = simulation
            .cube
            .places
            .traps
            .iter()
            .map(|point| (point.x, point.y, point.z))
            .collect();

        simulation.regenerate_cube().unwrap();
        let changed = original
            .iter()
            .zip(&simulation.cube.places.traps)
            .any(|(old, new)| old.0 != new.x || old.1 != new.y || old.2 != new.z);
        assert!(changed);
        assert_eq!(simulation.cube.places.traps.len(), original.len());
        assert!(simulation
            .cube
            .places
            .traps
            .iter()
            .all(|point| simulation.cube.contains(point)));

        simulation.time_temperature.advance(1.0);
        assert!(simulation.time_temperature.current_time() > 0.0);
        simulation.reset_time_temperature();
        assert_eq!(simulation.time_temperature.current_time(), 0.0);
        assert!((simulation.time_temperature.current_temperature() - 293.15).abs() < 1.0e-12);
    }
}
