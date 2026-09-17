// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Gillespie-style event selection for the standard trap-and-hole model.
//!
//! At each profile state, the engine constructs every enabled microscopic
//! pathway, samples an exponential waiting time for each independent 
//! process, and advances to the shortest lifetime. The temperature profile's
//! maximum step competes as an inert [`common::charge_transfer::Event::None`]
//! candidate so rates are
//! rebuilt whenever temperature changes by one kelvin or a control-point
//! boundary is reached.

use common::charge_transfer::{ElectronicState, Event, Candidate, TimedCandidate, RecordedEvent}; 
use common::place_ids::{PlaceAvailability, PlaceId};
use common::trap_hole_band_tail::{TrapParameterLayout, TrapParameters};
use common::crystal::Cube;

use common::rate_equation_selection::Transitions;
use common::rate_equations::{ground_excited_state_weights};
use common::time_temperature::TimeTemperature;
use common::trap_hole_band_tail::ElectronPlaces;

use common::numeric::{Float, TimeFloat};
use io::inputs::SimulationInputs;
use io::outputs::append_monte_carlo_experiment_batch_to_file;
use rand::Rng;
use std::path::Path;

use std::time::Instant;


/// Calculate thermal ground/excited occupation probabilities for one trap family.
///
/// The weights follow the two-state balance between Arrhenius excitation and
/// the configured excited-to-ground relaxation frequency.
fn state_weights(
    parameters: &TrapParameters,
    temperature: Float,
) -> Result<(TimeFloat, TimeFloat), String> {
    ground_excited_state_weights(
        &parameters.excited_energy_gap,
        &parameters.s_frequency_e,
        &parameters.s_frequency_g,
        &temperature,
    )
    .ok_or_else(|| "could not calculate ground/excited-state weights".to_string())
}

/// Construct the one system-wide irradiation-filling candidate.
///
/// Dose rate is converted from the configured denominator to dose per second.
/// The total hazard scales with the number of empty traps; a concrete target
/// is deliberately deferred until the aggregate event wins.
fn aggregate_filling_candidate(
    occupied_population: usize,
    total_population: usize,
    inputs: &SimulationInputs,
    transitions: &Transitions,
) -> Result<Candidate, String> {
    let characteristic_dose = *inputs
        .filling
        .d0
        .first()
        .ok_or_else(|| "filling d0 requires at least one value".to_string())?;
    let configured_dose_rate = *inputs
        .filling
        .d_dot
        .first()
        .ok_or_else(|| "filling d_dot requires at least one value".to_string())?;
    let seconds_per_dose_rate_unit =
        common::constants::time::unit_multiplier(inputs.filling.dd_unit)
            .ok_or_else(|| format!("unknown filling dose-rate unit: {}", inputs.filling.dd_unit))?
            .get_float_precision();
    let dose_rate = configured_dose_rate / seconds_per_dose_rate_unit;

    Candidate::filling_candidate(
        transitions.get_filling_transitions(),
        &characteristic_dose,
        &dose_rate,
        occupied_population,
        total_population,
    )
}

/// Sample all enabled events and update the current minimum lifetime.
///
/// Localised tunnelling candidates are created for every occupied-source/
/// available-destination pair using boundary-aware distances. Thermal release
/// contributes two candidates per occupied trap, and filling contributes one
/// population-level candidate.
/// Trial [Candidate] rates are calculated from [common::charge_transfer::Candidate]
/// These are converted to [TimedCandidate] and only retained if a shorter time is generated.
fn build_candidates(
    places: &ElectronPlaces,
    trap_places: &PlaceAvailability,
    hole_places: &PlaceAvailability,
    trap_parameters: &TrapParameterLayout,
    temperature: Float,
    cube: &Cube,
    inputs: &SimulationInputs,
    transitions: &Transitions,
    rng: &mut impl Rng,
    shortest: &mut TimedCandidate, 
) -> Result<(), String> {



    
    for &source in trap_places.available() {

        let parameters = trap_parameters.get(source);
        let (ground_weight, excited_weight) = state_weights(parameters, temperature)?;

        if transitions.get_delocalised(){
            let (ground, excited) = Candidate::delocalised_candidates(
                                                                      transitions.get_delocaised_transitions(), 
                                                                      &parameters, source, temperature, 
                                                                      ground_weight, excited_weight)?;
            shortest.find_shortest(ground, rng)?;
            shortest.find_shortest(excited, rng)?;
                                                
        }
       
        if transitions.get_localised_recombination(){
            for &hole in hole_places.available() {
                let distance = cube.distance(
                    &places.traps()[source.index()],
                    &places.holes()[hole.index()],
                );
                let (ground, excited) = Candidate::localised_recombination_candidates(
                                                                          transitions.get_locaised_recomb_transitions(), 
                                                                          &parameters, source, hole, temperature, 
                                                                          distance, ground_weight, excited_weight)?;
                shortest.find_shortest(ground, rng)?;
                shortest.find_shortest(excited, rng)?;
     
           } 
        }
        if transitions.get_localised_retrapping(){
            for &destination in trap_places.unavailable() {
                if destination == source{
                    continue
                }
                let distance = cube.distance(
                    &places.traps()[source.index()],
                    &places.traps()[destination.index()],
                );

                let (ground, excited) = Candidate::localised_retrapping_candidates(
                                                                          transitions.get_locaised_retrap_transitions(), 
                                                                          &parameters, source, destination, temperature, 
                                                                          distance, ground_weight, excited_weight)?;
                shortest.find_shortest(ground, rng)?;
                shortest.find_shortest(excited, rng)?;

            }
        }
    }
  
    if transitions.get_filling(){
        let fill = aggregate_filling_candidate(
            trap_places.available_count(),
            trap_places.total(),
            inputs,
            transitions,
        )?;
        shortest.find_shortest(fill, rng)?;
           
    }  

    Ok(())
}

/// Apply the selected event to trap and hole occupancy partitions.
///
/// Direct tunnelling has a destination already attached. A thermal-release
/// event first performs a second distance-weighted competition between
/// conduction-band recombination and retrapping destinations. The returned
/// event is the fully resolved transition stored in output.
fn apply_event(
    event: Event,
    places: &ElectronPlaces,
    trap_places: &mut PlaceAvailability,
    hole_places: &mut PlaceAvailability,
    trap_parameters: &TrapParameterLayout,
    cube: &Cube,
    transitions: &Transitions,
    rng: &mut impl Rng,
) -> Result<Event, String> {
    match event {
        Event::LocalisedRecombination { source, hole, .. }
        | Event::DelocalisedRecombination { source, hole, .. } => {
            if !trap_places.make_unavailable(source) {
                return Err(format!("recombination source {source:?} was not occupied"));
            }
            if !hole_places.make_unavailable(hole) {
                return Err(format!("recombination hole {hole:?} was not available"));
            }
        }
        Event::LocalisedRetrapping { source, destination, .. }
        | Event::DelocalisedRetrapping { source, destination, .. } => {
            if !trap_places.make_unavailable(source) {
                return Err(format!("retrapping source {source:?} was not occupied"));
            }
            if !trap_places.make_available(destination) {
                return Err(format!(
                    "retrapping destination {destination:?} was already occupied"
                ));
            }
        }
        Event::Delocalised { source, state } => {
            let outcome = choose_delocalised_outcome(
                source,
                places,
                trap_places,
                hole_places,
                trap_parameters.get(source),
                cube,
                transitions,
                state,
                rng
            )?;

            if !trap_places.make_unavailable(source) {
                return Err(format!("delocalised source {source:?} was not occupied"));
            }

            match outcome {
                TimedCandidate {
                    event: selected_event @ Event::DelocalisedRecombination { hole, .. },
                    ..
                } => {
                    if !hole_places.make_unavailable(hole) {
                        return Err(format!(
                            "delocalised recombination hole {hole:?} was not available"
                        ));
                    }
                    return Ok(selected_event);
                }
                TimedCandidate {
                    event: selected_event @ Event::DelocalisedRetrapping { destination, .. },
                    ..
                } => {
                    if !trap_places.make_available(destination) {
                        return Err(format!(
                            "delocalised retrapping destination {destination:?} was occupied"
                        ));
                    }
                    return Ok(selected_event);
                }
                _ => {
                    return Err(format!(
                            "An event that is not Delocalised Recombination or Retrapping \\
                            has been returned when a delocalised transition has been selected.
                            \\This error should never occur. Panic -- A LOT!"
                    ));
                }
            }
        }
        Event::Filling { .. } => {
            let filling = choose_filling_outcome(trap_places, hole_places, transitions, rng)?;
            return Ok(filling)
        }
        Event::None => return Ok(Event::None),
    }

    Ok(event)
}

/// Choose where a conduction-band electron is captured.
///
/// Each active hole and, when enabled, each empty trap receives a sampled
/// time proportional to `exp((r / mu)^2)`. The smallest time favours nearby
/// centres. The current reciprocal-rate convention multiplies trap waiting
/// times by `retrap_ratio`; larger positive ratios therefore make retrapping
/// slower relative to the unit recombination factor.
fn choose_delocalised_outcome(
    source: PlaceId,
    places: &ElectronPlaces,
    trap_places: &PlaceAvailability,
    hole_places: &PlaceAvailability,
    parameters: &TrapParameters,
    cube: &Cube,
    transitions: &Transitions,
    state: ElectronicState,
    rng: &mut impl Rng,
) -> Result<TimedCandidate, String> {
    let mu = parameters.delocalised_mu;
    let recombination_prefactor = 1.0;

    if !mu.is_finite() {
        return Err(format!(
            "delocalised mu must be finite and greater than zero, got {mu}"
        ));
    } else if mu <= 0.0 {
        if transitions.get_conduction_band_retrapping() && 
           trap_places.unavailable_count()> 0 && 
           rng.gen_bool(parameters.retrap_ratio/(1.0+parameters.retrap_ratio))
        {
            let trap = {
                let trap_dest = trap_places.unavailable();
                if trap_dest.is_empty() {
                    return Err("filling selected when no empty traps remain".to_string());
                }
                trap_dest[rng.gen_range(0..trap_dest.len())]
            };
        return Ok( TimedCandidate { 
            event: 
                Event::DelocalisedRetrapping { source, destination: trap, state }, 
            time: 0.0 });
        } else {
            let hole_destination = {
                let empty_holes = hole_places.available();
                if empty_holes.is_empty() {
                    return Err("No holes to put electron in".to_string());
                }
                empty_holes[rng.gen_range(0..empty_holes.len())]
            };
            return Ok( TimedCandidate { 
                event: 
                    Event::DelocalisedRecombination { source, hole: hole_destination, state }, 
                time: 0.0 });
            
        } 
    } else {
   
        let source_position = &places.traps()[source.index()];
    
        let mut current_shortest = TimedCandidate::new_negative_time();

        for &hole in hole_places.available() {
            let distance = cube.distance(source_position, &places.holes()[hole.index()]);
            current_shortest.find_smallest_candidate(TimedCandidate::delocalised_recombination(
                recombination_prefactor,
                mu,
                distance,
                source,
                hole,
                state,
                rng,
                )?
            )?;
        }

        if transitions.get_conduction_band_retrapping() && parameters.retrap_ratio > 0.0 {
            let retrapping_prefactor = recombination_prefactor*parameters.retrap_ratio;
            for &destination in trap_places.unavailable() {
                let distance = cube.distance(source_position, &places.traps()[destination.index()]);
                
                current_shortest.find_smallest_candidate(TimedCandidate::delocalised_retrapping(
                    retrapping_prefactor,
                    mu,
                    distance,
                    source,
                    destination,
                    state,
                    rng,
                    )?
                )?;
            }
        }
        Ok(current_shortest)
    }
    
}

/// Resolve an aggregate irradiation event into population changes.
///
/// Irradiation first activates one previously inactive hole site. Normally it
/// also occupies a uniformly selected empty trap. When filling-time
/// recombination is enabled, a fixed 0.5 branch instead consumes an active
/// hole; the returned [`Event::Filling`] records the selected identifiers.
pub fn choose_filling_outcome(
    trap_places: &mut PlaceAvailability,
    hole_places: &mut PlaceAvailability,
    transitions: &Transitions,
    rng: &mut impl Rng,
) -> Result<Event, String> {
    
    let hole = {
        let empty_holes = hole_places.unavailable();
        if empty_holes.is_empty() {
            return Err("filling selected when no available holes remain".to_string());
        }

        empty_holes[rng.gen_range(0..empty_holes.len())]
    };

    if !hole_places.make_available(hole) {
        return Err(format!("filling destination {hole:?} was occupied"));
    }

    if transitions.get_filling_retrapping() && rng.gen_bool(0.5){
        let hole_destination = {
            let empty_holes = hole_places.available();
            if empty_holes.is_empty() {
                return Err("No holes to put electron in".to_string());
            }
            empty_holes[rng.gen_range(0..empty_holes.len())]
        };
        if !hole_places.make_unavailable(hole_destination) {
            return Err(format!("hole destination {hole_destination:?} was occupied"));
        }
        return Ok(Event::Filling { trap:hole_destination, hole });
    } else { 
        let trap = {
            let trap_dest = trap_places.unavailable();
            if trap_dest.is_empty() {
                return Err("filling selected when no empty traps remain".to_string());
            }
            trap_dest[rng.gen_range(0..trap_dest.len())]
        };

        if !trap_places.make_available(trap) {
            return Err(format!("filling destination {trap:?} was occupied"));
        }
        return Ok(Event::Filling { trap, hole });
    }

}

/// Run one standard kinetic Monte Carlo trajectory to the profile endpoint.
///
/// Competing exponential lifetimes are resampled after every physical event
/// and temperature-profile boundary. The initial state and every resulting
/// state are written as [`RecordedEvent`] values in batches using
/// `results.capacity()` as the flush threshold.
pub fn run_standard(
    places: &ElectronPlaces,
    trap_places: &mut PlaceAvailability,
    hole_places: &mut PlaceAvailability,
    trap_parameters: &TrapParameterLayout,
    time_temperature: &mut TimeTemperature,
    cube: &Cube,
    inputs: &SimulationInputs,
    transitions: &Transitions,
    output_file: &Path,
    mut results: Vec<RecordedEvent>,
    rng: &mut impl Rng,
) -> Result<(), String> {
    results.push(RecordedEvent {
                    time: time_temperature.current_time(),
                    fill: trap_places.fill_ratio(),
                    temperature: time_temperature.current_temperature(),
                    event: Event::None,
                });
    while time_temperature.current_max_dt() != 0.0 {
        let temperature = time_temperature.current_temperature();

        // current_max_dt is signed because geological profiles run backwards.
        let signed_profile_dt = time_temperature.current_max_dt();
        let max_dt = signed_profile_dt.abs();
        let direction = signed_profile_dt.signum();

        let mut next_event = TimedCandidate { event: Event::None, time: max_dt,};
        // Contains localised, delocalised, and one aggregate filling event.
        build_candidates(
            places,
            trap_places,
            hole_places,
            trap_parameters,
            temperature,
            cube,
            inputs,
            transitions,
            rng,
            &mut next_event,
        )?;
    //     match next_event.event{
    //          Event::LocalisedRetrapping{source, destination, state} => {
    //              print!("Source: {}| destination: {}", source.index(), destination.index());
    //          },
    //          _ => {}
    //     };
        let signed_event_dt = direction * next_event.time;
        time_temperature.advance(signed_event_dt);
        let applied_event = apply_event(
                    next_event.event,
                    places,
                    trap_places,
                    hole_places,
                    trap_parameters,
                    cube,
                    transitions,
                    rng,
                )?;
                results.push(RecordedEvent {
                    time: time_temperature.current_time(),
                    fill: trap_places.fill_ratio(),
                    temperature: time_temperature.current_temperature(),
                    event: applied_event,
                });

        if results.len() == results.capacity() {
            append_monte_carlo_experiment_batch_to_file(output_file, &results)
                .map_err(|error| error.to_string())?;
            results.clear();
        }
    }

    if results.len() != 0 {
        append_monte_carlo_experiment_batch_to_file(output_file, &results)
            .map_err(|error| error.to_string())?;
        results.clear();
    }

    Ok(())
}



#[cfg(test)]
mod tests {
    use super::*;
    use common::trap_hole_band_tail::Coord;

    fn selected_transitions(
        localised_recombination: bool,
        localised_retrapping: bool,
        delocalised: bool,
        filling: bool,
    ) -> Transitions {
        Transitions::from_bool(
            localised_recombination,
            localised_recombination,
            delocalised,
            delocalised,
            filling,
            "first",
            localised_retrapping,
            localised_retrapping,
            false,
            false,
        )
        .unwrap()
    }

    fn parameters() -> TrapParameters {
        TrapParameters::new(
            0.1, 
            1.0e12, 
            1.0e12, 
            0.2, 
            0.1,
            1.0e12, 
            1.0e12, 
            1.0e12, 
            1.0e12, 
            1.0, 
            1.0, 
            0.1, 
            0.5,
        )
    }

    fn cb_retrapping_transitions() -> Transitions {
        Transitions::from_bool(
            false, false, true, false, false, "first", false, false, true, false,
        )
        .unwrap()
    }

    fn delocalised_test_state() -> (ElectronPlaces, PlaceAvailability, PlaceAvailability, Cube) {
        let places = ElectronPlaces::new_standard(
            vec![
                Coord::new(0.0, 0.0, 0.0).unwrap(),
                Coord::new(0.5, 0.0, 0.0).unwrap(),
            ],
            vec![Coord::new(0.25, 0.0, 0.0).unwrap()],
        )
        .unwrap();
        let mut trap_places = PlaceAvailability::new(2).unwrap();
        trap_places.make_available(PlaceId::new(0).unwrap());
        let mut hole_places = PlaceAvailability::new(1).unwrap();
        hole_places.make_available(PlaceId::new(0).unwrap());
        let cube = Cube::new(1.0, 1.0, 1.0, 2, 1, 0, false).unwrap();

        (places, trap_places, hole_places, cube)
    }

    #[test]
    fn candidate_builders_add_ground_and_excited_candidates() {
        let source = PlaceId::new(0).unwrap();
        let destination = PlaceId::new(1).unwrap();
        let hole = PlaceId::new(0).unwrap();
        let transitions = selected_transitions(true, true, true, false);
        let parameters = parameters();
        let temperature = 300.0;
        let (ground_weight, excited_weight) =
            state_weights(&parameters, temperature).unwrap();

        let delocalised = Candidate::delocalised_candidates(
            transitions.get_delocaised_transitions(),
            &parameters,
            source,
            temperature,
            ground_weight,
            excited_weight,
        )
        .unwrap();
        let recombination = Candidate::localised_recombination_candidates(
            transitions.get_locaised_recomb_transitions(),
            &parameters,
            source,
            hole,
            temperature,
            0.5,
            ground_weight,
            excited_weight,
        )
        .unwrap();
        let retrapping = Candidate::localised_retrapping_candidates(
            transitions.get_locaised_retrap_transitions(),
            &parameters,
            source,
            destination,
            temperature,
            0.5,
            ground_weight,
            excited_weight,
        )
        .unwrap();

        let candidates = [
            delocalised.0,
            delocalised.1,
            recombination.0,
            recombination.1,
            retrapping.0,
            retrapping.1,
        ];

        assert_eq!(candidates.len(), 6);
        assert!(candidates.iter().all(|candidate| candidate.rate > 0.0));
    }

    #[test]
    fn filling_is_one_aggregate_candidate_with_one_sampled_time() {
        let mut inputs = SimulationInputs::default();
        inputs.filling.d0 = vec![4.0];
        inputs.filling.d_dot = vec![2.0];
        inputs.filling.dd_unit = common::constants::time::TimeUnit::Minute;
        let transitions = selected_transitions(false, false, false, true);
        let candidate =
            aggregate_filling_candidate(1, 4, &inputs, &transitions).unwrap();

        assert_eq!(
            candidate.event,
            Event::Filling {
                trap: PlaceId::new(0).unwrap(),
                hole: PlaceId::new(0).unwrap(),
            }
        );
        assert_eq!(candidate.rate, 0.025);

        let mut rng = common::random::get_std_rng_for_rep(0);
        let timed = TimedCandidate::rate_to_lifetime(candidate, &mut rng).unwrap();
        assert_eq!(timed.event, candidate.event);
        assert!(timed.time.is_finite() && timed.time > 0.0);
    }

    #[test]
    fn fully_occupied_traps_produce_no_finite_filling_event() {
        let inputs = SimulationInputs::default();
        let transitions = selected_transitions(false, false, false, true);
        let candidate =
            aggregate_filling_candidate(4, 4, &inputs, &transitions).unwrap();
        assert_eq!(candidate.rate, 0.0);

        let mut rng = common::random::get_std_rng_for_rep(0);
        let timed = TimedCandidate::rate_to_lifetime(candidate, &mut rng).unwrap();
        assert!(timed.time.is_infinite());
    }

    #[test]
    fn delocalised_outcome_uses_retrap_ratio_and_updates_selected_destination() {
        let source = PlaceId::new(0).unwrap();
        let destination = PlaceId::new(1).unwrap();
        let hole = PlaceId::new(0).unwrap();
        let transitions = cb_retrapping_transitions();

        let (places, mut trap_places, mut hole_places, cube) = delocalised_test_state();
        let mut recombination_parameters = parameters();
        recombination_parameters.delocalised_mu = 1.0;
        recombination_parameters.retrap_ratio = 0.0;
        let recombination_layout = TrapParameterLayout::uniform(recombination_parameters);
        let mut rng = common::random::get_std_rng_for_rep(0);

        let event = apply_event(
            Event::Delocalised {
                source,
                state: ElectronicState::Ground,
            },
            &places,
            &mut trap_places,
            &mut hole_places,
            &recombination_layout,
            &cube,
            &transitions,
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            event,
            Event::DelocalisedRecombination {
                source,
                hole,
                state: ElectronicState::Ground,
            }
        );
        assert_eq!(trap_places.available_count(), 0);
        assert_eq!(hole_places.available_count(), 0);

        let (places, mut trap_places, mut hole_places, cube) = delocalised_test_state();
        let mut retrapping_parameters = parameters();
        retrapping_parameters.delocalised_mu = 1.0;
        retrapping_parameters.retrap_ratio = 1.0;
        let retrapping_layout = TrapParameterLayout::uniform(retrapping_parameters);

        let event = apply_event(
            Event::Delocalised {
                source,
                state: ElectronicState::Excited,
            },
            &places,
            &mut trap_places,
            &mut hole_places,
            &retrapping_layout,
            &cube,
            &transitions,
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            event,
            Event::DelocalisedRecombination  {
                source,
                hole,
                state: ElectronicState::Excited,
            }
        );
        assert_eq!(trap_places.available_count(), 0);
        assert_eq!(hole_places.available_count(), 0);
    }
}
