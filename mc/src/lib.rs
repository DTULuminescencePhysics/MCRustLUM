// //! Kinetic Monte Carlo model setup and event-time generation.

// pub mod calculate_times;
// pub mod system_setup;

// pub use calculate_times::{
//     CalculationError, CompletedMonteCarloEvent, CompletedMonteCarloStep, MonteCarloEvent,
//     MonteCarloStep,
// };
// pub use system_setup::MonteCarloTrial;

// use common::crystal::Cube;
// use common::rate_equation_selection::Transitions as RateEquationTransitions;
// use common::time_temperature::TimeTemperature;
// use system_setup::{FillingParameters, SetupError, TrapParameters};

// /// Reusable geometry and equation configuration for Monte Carlo trials.
// ///
// /// The crystal and time-temperature profile are retained here so multiple
// /// trials can reuse them. The profile is reset whenever a trial is prepared.
// /// Per-trial parameters, distances, and occupancy are held in
// /// [`MonteCarloTrial`]. A trial is optional because the reusable model can be
// /// created before its site parameters are available.
// #[derive(Debug, Clone)]
// pub struct MonteCarlo {
//     cube: Cube,
//     rate_equations: RateEquationTransitions,
//     time_temperature: TimeTemperature,
//     monte_carlo_trial: Option<MonteCarloTrial>,
// }

// impl MonteCarlo {
//     /// Create reusable Monte Carlo setup without allocating trial distances.
//     ///
//     /// `time_temperature` is owned by the system and reused for every trial.
//     pub fn new(
//         cube: Cube,
//         rate_equations: RateEquationTransitions,
//         time_temperature: TimeTemperature,
//     ) -> Self {
//         Self {
//             cube,
//             rate_equations,
//             time_temperature,
//             monte_carlo_trial: None,
//         }
//     }

//     /// Crystal geometry from which the next trial will be prepared.
//     pub fn cube(&self) -> &Cube {
//         &self.cube
//     }

//     /// Equation selection shared by every trial until it is replaced.
//     pub fn rate_equations(&self) -> &RateEquationTransitions {
//         &self.rate_equations
//     }

//     /// Time and temperature profile shared by successive trials.
//     pub fn time_temperature(&self) -> &TimeTemperature {
//         &self.time_temperature
//     }

//     /// Mutable profile access for advancing a Monte Carlo calculation.
//     pub fn time_temperature_mut(&mut self) -> &mut TimeTemperature {
//         &mut self.time_temperature
//     }

//     /// Current trial, if its parameters and distances have been prepared.
//     pub fn monte_carlo_trial(&self) -> Option<&MonteCarloTrial> {
//         self.monte_carlo_trial.as_ref()
//     }

//     /// Mutable access to the current trial's occupancy and event operations.
//     pub fn monte_carlo_trial_mut(&mut self) -> Option<&mut MonteCarloTrial> {
//         self.monte_carlo_trial.as_mut()
//     }

//     /// Build the trial by referencing the currently stored crystal.
//     ///
//     /// The crystal is borrowed only while distances are calculated. The
//     /// resulting trial does not retain a crystal reference or clone.
//     pub fn set_monte_carlo_trial(
//         &mut self,
//         trap_parameters: TrapParameters,
//         filling_parameters: FillingParameters,
//     ) -> Result<&MonteCarloTrial, SetupError> {
//         let trial = MonteCarloTrial::from_crystal(
//             &self.cube,
//             &self.rate_equations,
//             trap_parameters,
//             filling_parameters,
//         )?;
//         self.time_temperature.reset();
//         self.monte_carlo_trial = Some(trial);
//         Ok(self
//             .monte_carlo_trial
//             .as_ref()
//             .expect("the trial was assigned immediately above"))
//     }

//     /// Replace the crystal and discard trial data tied to the old geometry.
//     pub fn replace_cube(&mut self, cube: Cube) {
//         self.cube = cube;
//         self.monte_carlo_trial = None;
//     }

//     /// Replace equation selection and discard pathway-dependent distance data.
//     pub fn replace_rate_equations(&mut self, rate_equations: RateEquationTransitions) {
//         self.rate_equations = rate_equations;
//         self.monte_carlo_trial = None;
//     }

//     /// Start a new trial with a replacement crystal in one consistent update.
//     ///
//     /// Trial construction is completed before stored state is changed, so an
//     /// error leaves the previous crystal, profile position, and trial untouched.
//     pub fn start_new_trial(
//         &mut self,
//         cube: Cube,
//         trap_parameters: TrapParameters,
//         filling_parameters: FillingParameters,
//     ) -> Result<&MonteCarloTrial, SetupError> {
//         let trial = MonteCarloTrial::from_crystal(
//             &cube,
//             &self.rate_equations,
//             trap_parameters,
//             filling_parameters,
//         )?;
//         self.time_temperature.reset();
//         self.cube = cube;
//         self.monte_carlo_trial = Some(trial);
//         Ok(self
//             .monte_carlo_trial
//             .as_ref()
//             .expect("the trial was assigned immediately above"))
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use common::rate_equation_selection::{
//         DelocalisedRateEquation, DelocalisedRateEquationType, FillingRateEquation,
//         LocalisedRateEquation, TransitionsTypes,
//     };
//     use system_setup::{
//         ConductionBandParameters, ConductionBandPathwayParameters, DelocalisedParameters,
//         LocalisedParameters, StateParameters, TrapParameter,
//     };

//     fn state(value: f32) -> StateParameters<TrapParameter<f32>> {
//         StateParameters {
//             ground: TrapParameter::Shared(value),
//             excited: TrapParameter::Shared(value),
//         }
//     }

//     fn trap_parameters(trap_count: usize) -> TrapParameters {
//         TrapParameters {
//             delocalised: DelocalisedParameters {
//                 energy_to_conduction_band: StateParameters {
//                     ground: TrapParameter::PerTrap(vec![0.3; trap_count]),
//                     excited: TrapParameter::Shared(0.2),
//                 },
//                 frequency: state(1.0e12),
//             },
//             localised_recombination: LocalisedParameters {
//                 alpha: state(1.0),
//                 frequency: state(2.0),
//             },
//             localised_retrapping: LocalisedParameters {
//                 alpha: state(0.8),
//                 frequency: state(3.0),
//             },
//             conduction_band: ConductionBandParameters {
//                 recombination: ConductionBandPathwayParameters {
//                     constant: 2.0,
//                     distance_scale: 1.0,
//                 },
//                 retrapping: ConductionBandPathwayParameters {
//                     constant: 3.0,
//                     distance_scale: 1.5,
//                 },
//             },
//             excited_energy_gap: TrapParameter::Shared(0.045),
//         }
//     }

//     fn filling_parameters() -> FillingParameters {
//         FillingParameters {
//             characteristic_dose: 4.0,
//             dose_rate: 2.0,
//         }
//     }

//     fn rate_equations() -> RateEquationTransitions {
//         RateEquationTransitions::NoCbFillRetrapping {
//             transitions: TransitionsTypes {
//                 delocalised: DelocalisedRateEquation::Both {
//                     re: DelocalisedRateEquationType::FirstOrder,
//                 },
//                 localised_recomb: LocalisedRateEquation::Both,
//                 localised_retrap: LocalisedRateEquation::None,
//                 filling: FillingRateEquation::Basic,
//             },
//         }
//     }

//     fn time_temperature() -> TimeTemperature {
//         TimeTemperature::new(vec![0.0, 10.0, 20.0], vec![300.0, 400.0, 500.0], "s").unwrap()
//     }

//     #[test]
//     fn model_creates_trial_from_a_borrowed_crystal() {
//         let cube = Cube::new_random(10.0, 10.0, 10.0, 3, 2, 0, false);
//         let mut model = MonteCarlo::new(cube, rate_equations(), time_temperature());

//         assert!(model.monte_carlo_trial().is_none());
//         let trial = model
//             .set_monte_carlo_trial(trap_parameters(3), filling_parameters())
//             .unwrap();

//         assert_eq!(trial.trap_count(), 3);
//         assert_eq!(trial.hole_count(), 2);
//         assert_eq!(trial.trap_occupancy(), &[false, false, false]);
//         assert_eq!(trial.hole_availability(), &[true, true]);
//         assert_eq!(model.cube().places.traps.len(), 3);
//     }

//     #[test]
//     fn replacing_crystal_invalidates_stale_trial() {
//         let cube = Cube::new_random(10.0, 10.0, 10.0, 2, 1, 0, false);
//         let mut model = MonteCarlo::new(cube, rate_equations(), time_temperature());
//         model
//             .set_monte_carlo_trial(trap_parameters(2), filling_parameters())
//             .unwrap();

//         model.replace_cube(Cube::new_random(20.0, 20.0, 20.0, 4, 3, 0, true));

//         assert!(model.monte_carlo_trial().is_none());
//         assert_eq!(model.cube().places.traps.len(), 4);
//     }

//     #[test]
//     fn new_trial_updates_crystal_and_trial_together() {
//         let cube = Cube::new_random(10.0, 10.0, 10.0, 2, 1, 0, false);
//         let mut model = MonteCarlo::new(cube, rate_equations(), time_temperature());
//         model
//             .set_monte_carlo_trial(trap_parameters(2), filling_parameters())
//             .unwrap();

//         let replacement = Cube::new_random(15.0, 15.0, 15.0, 5, 4, 0, true);
//         let trial = model
//             .start_new_trial(replacement, trap_parameters(5), filling_parameters())
//             .unwrap();

//         assert_eq!(trial.trap_count(), 5);
//         assert_eq!(trial.hole_count(), 4);
//         assert_eq!(model.cube().places.traps.len(), 5);
//         assert!(model.cube().boundary.kind == common::crystal::BoundaryCondition::Periodic);
//     }

//     #[test]
//     fn preparing_each_trial_resets_the_shared_time_temperature_profile() {
//         let cube = Cube::new_random(10.0, 10.0, 10.0, 2, 1, 0, false);
//         let mut model = MonteCarlo::new(cube, rate_equations(), time_temperature());

//         model.time_temperature_mut().advance(5.0);
//         assert_eq!(model.time_temperature().current_time(), 5.0);
//         assert_eq!(model.time_temperature().current_temperature(), 350.0);

//         model
//             .set_monte_carlo_trial(trap_parameters(2), filling_parameters())
//             .unwrap();
//         assert_eq!(model.time_temperature().current_time(), 0.0);
//         assert_eq!(model.time_temperature().current_temperature(), 300.0);

//         model.time_temperature_mut().advance(10.0);
//         let replacement = Cube::new_random(15.0, 15.0, 15.0, 3, 2, 0, true);
//         model
//             .start_new_trial(replacement, trap_parameters(3), filling_parameters())
//             .unwrap();
//         assert_eq!(model.time_temperature().current_time(), 0.0);
//         assert_eq!(model.time_temperature().current_temperature(), 300.0);
//     }
// }
