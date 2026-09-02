// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only
use crate::experiment::{PlaceAvailability, PlaceId, TrapParameterLayout, TrapParameters};
use common::crystal::Cube;
use common::rate_equation_inputs::{
    DelocalisedTransitionInputs, FillingTransitionInputs, LocalisedTransitionInputs,
};
use common::rate_equation_selection::{
    DelocalisedRateEquation, FillingRateEquation, LocalisedRateEquation, Transitions,
    TransitionsTypes,
};
use common::rate_equations::{ground_excited_state_weights, retrapping_probability_by_r};
use common::time_temperature::TimeTemperature;
use common::trap_hole_band_tail::ElectronPlaces;

use common::numeric::{Float, TimeFloat};
use io::SimulationInputs;
use io::outputs::append_monte_carlo_experiment_batch_to_file;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElectronicState {
    Ground,
    Excited,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Event {
    LocalisedRecombination {
        source: PlaceId,
        hole: PlaceId,
        state: ElectronicState,
    },
    LocalisedRetrapping {
        source: PlaceId,
        destination: PlaceId,
        state: ElectronicState,
    },
    Delocalised {
        source: PlaceId,
        state: ElectronicState,
    },
    DelocalisedRecombination {
        source: PlaceId,
        hole: PlaceId,
        state: ElectronicState,
    },
    DelocalisedRetrapping {
        source: PlaceId,
        destination: PlaceId,
        state: ElectronicState,
    },
    Filling {
        trap: PlaceId,
        hole: PlaceId,
    },
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DelocalisedOutcome {
    Recombination { hole: PlaceId },
    Retrapping { destination: PlaceId },
}

#[derive(Debug, Clone, Copy)]
struct TimedDelocalisedOutcome {
    outcome: DelocalisedOutcome,
    time: TimeFloat,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    event: Event,
    rate: TimeFloat,
}

fn transition_types(transitions: &Transitions) -> &TransitionsTypes {
    match transitions {
        Transitions::NoCbFillRetrapping { transitions }
        | Transitions::FillRetrapping { transitions }
        | Transitions::CbRetrapping { transitions }
        | Transitions::FillCbRetrapping { transitions } => transitions,
    }
}

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

fn push_candidate(
    candidates: &mut Vec<Candidate>,
    event: Event,
    rate: Option<TimeFloat>,
) -> Result<(), String> {
    let rate = rate.ok_or_else(|| format!("could not calculate rate for {event:?}"))?;

    if !rate.is_finite() {
        return Err(format!("non-finite rate for {event:?}: {rate}"));
    }
    if rate < 0.0 {
        return Err(format!("negative rate for {event:?}: {rate}"));
    }
    if rate > 0.0 {
        candidates.push(Candidate { event, rate });
    }

    Ok(())
}

fn push_delocalised_candidates(
    candidates: &mut Vec<Candidate>,
    source: PlaceId,
    parameters: &TrapParameters,
    temperature: Float,
    transitions: &Transitions,
) -> Result<(), String> {
    if matches!(
        transition_types(transitions).delocalised,
        DelocalisedRateEquation::None
    ) {
        return Ok(());
    }

    let (ground_weight, excited_weight) = state_weights(parameters, temperature)?;
    let inputs = DelocalisedTransitionInputs {
        e_cb_ground: parameters.e_cb_ground,
        frequency_ground: parameters.de_frequency_ground,
        e_cb_excited: parameters.e_cb_excited,
        frequency_excited: parameters.de_frequency_excited,
        temperature,
        ground_weight,
        excited_weight,
    };
    let (ground_rate, excited_rate): (Option<TimeFloat>, Option<TimeFloat>) =
        transition_types(transitions).delocalised.calculate(&inputs);

    push_candidate(
        candidates,
        Event::Delocalised {
            source,
            state: ElectronicState::Ground,
        },
        ground_rate,
    )?;
    push_candidate(
        candidates,
        Event::Delocalised {
            source,
            state: ElectronicState::Excited,
        },
        excited_rate,
    )
}

fn push_recombination_candidates(
    candidates: &mut Vec<Candidate>,
    source: PlaceId,
    hole: PlaceId,
    distance: Float,
    parameters: &TrapParameters,
    temperature: Float,
    transitions: &Transitions,
) -> Result<(), String> {
    if matches!(
        transition_types(transitions).localised_recomb,
        LocalisedRateEquation::None
    ) {
        return Ok(());
    }

    let (ground_weight, excited_weight) = state_weights(parameters, temperature)?;
    let inputs = LocalisedTransitionInputs {
        alpha_ground: parameters.alpha_ground,
        frequency_ground: parameters.lo_frequency_ground,
        alpha_excited: parameters.alpha_excited,
        frequency_excited: parameters.lo_frequency_excited,
        ground_weight,
        excited_weight,
        distance,
    };
    let (ground_rate, excited_rate): (Option<TimeFloat>, Option<TimeFloat>) =
        transition_types(transitions)
            .localised_recomb
            .calculate(&inputs);

    push_candidate(
        candidates,
        Event::LocalisedRecombination {
            source,
            hole,
            state: ElectronicState::Ground,
        },
        ground_rate,
    )?;
    push_candidate(
        candidates,
        Event::LocalisedRecombination {
            source,
            hole,
            state: ElectronicState::Excited,
        },
        excited_rate,
    )
}

fn push_retrapping_candidates(
    candidates: &mut Vec<Candidate>,
    source: PlaceId,
    destination: PlaceId,
    distance: Float,
    parameters: &TrapParameters,
    temperature: Float,
    transitions: &Transitions,
) -> Result<(), String> {
    if matches!(
        transition_types(transitions).localised_retrap,
        LocalisedRateEquation::None
    ) {
        return Ok(());
    }

    let (ground_weight, excited_weight) = state_weights(parameters, temperature)?;
    let inputs = LocalisedTransitionInputs {
        alpha_ground: parameters.alpha_ground,
        frequency_ground: parameters.lo_frequency_ground,
        alpha_excited: parameters.alpha_excited,
        frequency_excited: parameters.lo_frequency_excited,
        ground_weight,
        excited_weight,
        distance,
    };
    let (ground_rate, excited_rate): (Option<TimeFloat>, Option<TimeFloat>) =
        transition_types(transitions)
            .localised_retrap
            .calculate(&inputs);

    push_candidate(
        candidates,
        Event::LocalisedRetrapping {
            source,
            destination,
            state: ElectronicState::Ground,
        },
        ground_rate,
    )?;
    push_candidate(
        candidates,
        Event::LocalisedRetrapping {
            source,
            destination,
            state: ElectronicState::Excited,
        },
        excited_rate,
    )
}

fn push_aggregate_filling_candidate(
    candidates: &mut Vec<Candidate>,
    occupied_population: usize,
    total_population: usize,
    inputs: &SimulationInputs,
    transitions: &Transitions,
) -> Result<(), String> {
    if matches!(
        transition_types(transitions).filling,
        FillingRateEquation::None
    ) {
        return Ok(());
    }

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
    let filling_inputs = FillingTransitionInputs {
        characteristic_dose,
        dose_rate,
        occupied_population: occupied_population as Float,
        total_population: total_population as Float,
    };
    let rate: Option<TimeFloat> = transition_types(transitions)
        .filling
        .calculate(&filling_inputs);

    push_candidate(
        candidates,
        Event::Filling {
            trap: PlaceId::new(0)?,
            hole: PlaceId::new(0)?,
        },
        rate,
    )
}

fn build_candidates(
    places: &ElectronPlaces,
    trap_places: &PlaceAvailability,
    hole_places: &PlaceAvailability,
    trap_parameters: &TrapParameterLayout,
    temperature: Float,
    cube: &Cube,
    inputs: &SimulationInputs,
    transitions: &Transitions,
) -> Result<Vec<Candidate>, String> {
    let mut candidates = Vec::new();

    for &source in trap_places.available() {
        let parameters = trap_parameters.get(source);

        // Delocalised rates: once per occupied source trap.
        push_delocalised_candidates(
            &mut candidates,
            source,
            parameters,
            temperature,
            transitions,
        )?;

        // Recombination: once per occupied-source/active-hole pair.
        for &hole in hole_places.available() {
            let distance = cube.distance(
                &places.traps()[source.index()],
                &places.holes()[hole.index()],
            );

            push_recombination_candidates(
                &mut candidates,
                source,
                hole,
                distance,
                parameters,
                temperature,
                transitions,
            )?;
        }

        // Retrapping: once per occupied-source/empty-destination pair.
        for &destination in trap_places.unavailable() {
            let distance = cube.distance(
                &places.traps()[source.index()],
                &places.traps()[destination.index()],
            );

            push_retrapping_candidates(
                &mut candidates,
                source,
                destination,
                distance,
                parameters,
                temperature,
                transitions,
            )?;
        }
    }

    // The filling rate already contains the number of empty traps.
    push_aggregate_filling_candidate(
        &mut candidates,
        trap_places.available_count(),
        places.traps().len(),
        inputs,
        transitions,
    )?;

    Ok(candidates)
}

#[derive(Debug, Clone, Copy)]
struct TimedCandidate {
    event: Event,
    time: TimeFloat,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RecordedEvent {
    pub time: TimeFloat,
    pub fill: Float,
    pub temperature: Float,
    pub event: Event,
}

fn lifetime<R: Rng + ?Sized>(rate: TimeFloat, rng: &mut R) -> Result<TimeFloat, String> {
    if !rate.is_finite() {
        return Err(format!("non-finite transition rate: {rate}"));
    }

    if rate < 0.0 {
        return Err(format!("negative transition rate: {rate}"));
    }

    if rate == 0.0 {
        return Ok(TimeFloat::INFINITY);
    }

    let u: TimeFloat = rng.sample(rand::distributions::Open01);
    Ok(-u.ln() / rate)
}

fn cb_retrapping_enabled(transitions: &Transitions) -> bool {
    matches!(
        transitions,
        Transitions::CbRetrapping { .. } | Transitions::FillCbRetrapping { .. }
    )
}

fn reciprocal_rate_lifetime<R: Rng + ?Sized>(
    reciprocal_rate: TimeFloat,
    rng: &mut R,
) -> Result<TimeFloat, String> {
    if reciprocal_rate.is_nan() || reciprocal_rate < 0.0 {
        return Err(format!(
            "invalid delocalised destination reciprocal rate: {reciprocal_rate}"
        ));
    }
    if reciprocal_rate.is_infinite() {
        return Ok(TimeFloat::INFINITY);
    }

    let u: TimeFloat = rng.sample(rand::distributions::Open01);
    Ok(-u.ln() * reciprocal_rate)
}

fn push_delocalised_outcome<R: Rng + ?Sized>(
    candidates: &mut Vec<TimedDelocalisedOutcome>,
    outcome: DelocalisedOutcome,
    prefactor: Float,
    mu: Float,
    distance: Float,
    rng: &mut R,
) -> Result<(), String> {
    let reciprocal_rate: Option<TimeFloat> =
        retrapping_probability_by_r(&prefactor, &mu, &distance);
    let reciprocal_rate = reciprocal_rate.ok_or_else(|| {
        format!("could not calculate delocalised destination rate for {outcome:?}")
    })?;
    let time = reciprocal_rate_lifetime(reciprocal_rate, rng)?;

    if time.is_finite() {
        candidates.push(TimedDelocalisedOutcome { outcome, time });
    }

    Ok(())
}

fn choose_delocalised_outcome(
    source: PlaceId,
    places: &ElectronPlaces,
    trap_places: &PlaceAvailability,
    hole_places: &PlaceAvailability,
    parameters: &TrapParameters,
    cube: &Cube,
    transitions: &Transitions,
) -> Result<Option<DelocalisedOutcome>, String> {
    let mu = parameters.delocalised_mu;
    let recombination_prefactor = parameters.retrap_ratio;

    if !mu.is_finite() || mu <= 0.0 {
        return Err(format!(
            "delocalised mu must be finite and greater than zero, got {mu}"
        ));
    }
    if !recombination_prefactor.is_finite() || !(0.0..=1.0).contains(&recombination_prefactor) {
        return Err(format!(
            "delocalised retrap_ratio must be between zero and one, got {recombination_prefactor}"
        ));
    }

    let retrapping_prefactor = 1.0 - recombination_prefactor;
    let source_position = &places.traps()[source.index()];
    let mut rng = common::random::rng();
    let mut candidates = Vec::with_capacity(
        hole_places.available().len()
            + if cb_retrapping_enabled(transitions) {
                trap_places.unavailable().len()
            } else {
                0
            },
    );

    for &hole in hole_places.available() {
        let distance = cube.distance(source_position, &places.holes()[hole.index()]);
        push_delocalised_outcome(
            &mut candidates,
            DelocalisedOutcome::Recombination { hole },
            recombination_prefactor,
            mu,
            distance,
            &mut *rng,
        )?;
    }

    if cb_retrapping_enabled(transitions) {
        for &destination in trap_places.unavailable() {
            let distance = cube.distance(source_position, &places.traps()[destination.index()]);
            push_delocalised_outcome(
                &mut candidates,
                DelocalisedOutcome::Retrapping { destination },
                retrapping_prefactor,
                mu,
                distance,
                &mut *rng,
            )?;
        }
    }

    let Some(minimum_time) = candidates
        .iter()
        .map(|candidate| candidate.time)
        .min_by(TimeFloat::total_cmp)
    else {
        return Ok(None);
    };

    // Zero prefactors can create exact ties. Pick uniformly among tied
    // destinations rather than depending on storage order.
    let tied_candidates: Vec<_> = candidates
        .iter()
        .filter(|candidate| candidate.time.total_cmp(&minimum_time).is_eq())
        .collect();
    let selected = tied_candidates[rng.gen_range(0..tied_candidates.len())];

    Ok(Some(selected.outcome))
}

fn calculate_candidate_times(
    candidates: &[Candidate],
    max_dt: TimeFloat,
) -> Result<Vec<TimedCandidate>, String> {
    let mut rng = common::random::rng();
    let mut timed = Vec::with_capacity(candidates.len());

    for candidate in candidates {
        timed.push(TimedCandidate {
            event: candidate.event,
            time: lifetime(candidate.rate, &mut *rng)?,
        });
    }
    timed.push(TimedCandidate {
        event: Event::None,
        time: max_dt,
    });
    Ok(timed)
}

fn earliest_candidate(candidates: &[TimedCandidate]) -> Option<TimedCandidate> {
    candidates
        .iter()
        .copied()
        .min_by(|a, b| a.time.total_cmp(&b.time))
}

pub(crate) fn run_standard(
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
) -> Result<(), String> {
    while time_temperature.current_max_dt() != 0.0 {
        let temperature = time_temperature.current_temperature();

        // current_max_dt is signed because geological profiles run backwards.
        let signed_profile_dt = time_temperature.current_max_dt();
        let max_dt = signed_profile_dt.abs();
        let direction = signed_profile_dt.signum();

        // Contains localised, delocalised, and one aggregate filling event.
        let rate_candidates = build_candidates(
            places,
            trap_places,
            hole_places,
            trap_parameters,
            temperature,
            cube,
            inputs,
            transitions,
        )?;

        let timed_candidates = calculate_candidate_times(&rate_candidates, max_dt)?;

        let earliest = earliest_candidate(&timed_candidates);

        match earliest {
            Some(TimedCandidate { event, time }) => {
                let signed_event_dt = direction * time;
                time_temperature.advance(signed_event_dt);
                let applied_event = apply_event(
                    event,
                    places,
                    trap_places,
                    hole_places,
                    trap_parameters,
                    cube,
                    transitions,
                )?;
                results.push(RecordedEvent {
                    time: time_temperature.current_time(),
                    fill: trap_places.fill_ratio(),
                    temperature: time_temperature.current_temperature(),
                    event: applied_event,
                });
            }
            None => {
                time_temperature.advance(time_temperature.current_max_dt());
                results.push(RecordedEvent {
                    time: time_temperature.current_time(),
                    fill: trap_places.fill_ratio(),
                    temperature: time_temperature.current_temperature(),
                    event: Event::None,
                });
            }
        }

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

fn apply_event(
    event: Event,
    places: &ElectronPlaces,
    trap_places: &mut PlaceAvailability,
    hole_places: &mut PlaceAvailability,
    trap_parameters: &TrapParameterLayout,
    cube: &Cube,
    transitions: &Transitions,
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
        Event::LocalisedRetrapping {
            source,
            destination,
            ..
        }
        | Event::DelocalisedRetrapping {
            source,
            destination,
            ..
        } => {
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
            )?;

            if !trap_places.make_unavailable(source) {
                return Err(format!("delocalised source {source:?} was not occupied"));
            }

            match outcome {
                Some(DelocalisedOutcome::Recombination { hole }) => {
                    if !hole_places.make_unavailable(hole) {
                        return Err(format!(
                            "delocalised recombination hole {hole:?} was not available"
                        ));
                    }
                    return Ok(Event::DelocalisedRecombination {
                        source,
                        hole,
                        state,
                    });
                }
                Some(DelocalisedOutcome::Retrapping { destination }) => {
                    if !trap_places.make_available(destination) {
                        return Err(format!(
                            "delocalised retrapping destination {destination:?} was occupied"
                        ));
                    }
                    return Ok(Event::DelocalisedRetrapping {
                        source,
                        destination,
                        state,
                    });
                }
                None => {}
            }
        }
        Event::Filling { trap: _, hole: _ } => {
            let trap = {
                let empty_traps = trap_places.unavailable();
                if empty_traps.is_empty() {
                    return Err("filling selected when no empty traps remain".to_string());
                }

                let mut rng = common::random::rng();
                empty_traps[rng.gen_range(0..empty_traps.len())]
            };

            if !trap_places.make_available(trap) {
                return Err(format!("filling destination {trap:?} was occupied"));
            }

            let hole = {
                let empty_holes = hole_places.unavailable();
                if empty_holes.is_empty() {
                    return Err("filling selected when no available holes remain".to_string());
                }

                let mut rng = common::random::rng();
                empty_holes[rng.gen_range(0..empty_holes.len())]
            };

            if !hole_places.make_available(hole) {
                return Err(format!("filling destination {hole:?} was occupied"));
            }
            return Ok(Event::Filling { trap, hole });
        }
        Event::None => return Ok(Event::None),
    }

    Ok(event)
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
            0.1, 1.0e12, 1.0e12, 0.2, 0.1, 1.0e12, 1.0e12, 1.0e12, 1.0e12, 1.0, 1.0, 0.1, 0.5,
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
    fn pathway_helpers_add_ground_and_excited_candidates() {
        let source = PlaceId::new(0).unwrap();
        let destination = PlaceId::new(1).unwrap();
        let hole = PlaceId::new(0).unwrap();
        let transitions = selected_transitions(true, true, true, false);
        let parameters = parameters();
        let mut candidates = Vec::new();

        push_delocalised_candidates(&mut candidates, source, &parameters, 300.0, &transitions)
            .unwrap();
        push_recombination_candidates(
            &mut candidates,
            source,
            hole,
            0.5,
            &parameters,
            300.0,
            &transitions,
        )
        .unwrap();
        push_retrapping_candidates(
            &mut candidates,
            source,
            destination,
            0.5,
            &parameters,
            300.0,
            &transitions,
        )
        .unwrap();

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
        let mut candidates = Vec::new();

        push_aggregate_filling_candidate(&mut candidates, 1, 4, &inputs, &transitions).unwrap();

        assert_eq!(candidates.len(), 1);
        assert_eq!(
            candidates[0].event,
            Event::Filling {
                trap: PlaceId::new(0).unwrap(),
                hole: PlaceId::new(0).unwrap(),
            }
        );
        assert_eq!(candidates[0].rate, 0.025);

        let timed = calculate_candidate_times(&candidates, 1.0).unwrap();
        assert_eq!(timed.len(), 2);
        assert_eq!(timed[0].event, candidates[0].event);
        assert!(timed[0].time.is_finite() && timed[0].time > 0.0);
        assert_eq!(timed[1].event, Event::None);
        assert_eq!(timed[1].time, 1.0);
    }

    #[test]
    fn filling_adds_no_candidate_when_every_trap_is_occupied() {
        let inputs = SimulationInputs::default();
        let transitions = selected_transitions(false, false, false, true);
        let mut candidates = Vec::new();

        push_aggregate_filling_candidate(&mut candidates, 4, 4, &inputs, &transitions).unwrap();

        assert!(candidates.is_empty());
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
        )
        .unwrap();

        assert_eq!(
            event,
            Event::DelocalisedRetrapping {
                source,
                destination,
                state: ElectronicState::Excited,
            }
        );
        assert_eq!(trap_places.available_count(), 1);
        assert!(trap_places.is_available(destination));
        assert_eq!(hole_places.available_count(), 1);
    }
}
