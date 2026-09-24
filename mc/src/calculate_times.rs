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

use common::charge_transfer::{ElectronicState, Event, Candidate, ReciprocalCandidate, TimedCandidate, RecordedEvent};
use common::charge_transfer::{delocalised_candidates, localised_recombination_candidates,localised_retrapping_candidates, 
    filling_candidate, gaussian_kernel, find_shortest_from_summed_reciprocal_rates}; 
use common::place_ids::{PlaceAvailability, PlaceId};
use common::trap_hole_band_tail::{TrapParameterLayout, ReTrapParameterLayout};
use common::crystal::Cube;
use common::rate_equation_inputs::FillingTransitionInputs;
use common::rate_equation_selection::Transitions;
use common::time_temperature::TimeTemperature;
use common::trap_hole_band_tail::{ElectronPlaces, Coord};

use common::numeric::Float;
use io::outputs::append_monte_carlo_experiment_batch_to_file;
use rand::{Rng, RngExt};
use std::path::Path;


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
    filling_inputs: &mut FillingTransitionInputs,
    transitions: &Transitions,
    rng: &mut impl Rng,
    shortest: &mut TimedCandidate, 
) -> Result<(), String> {

    
    if trap_places.available_count() > 0 {
        let (ground_weights, 
            excited_weights) = trap_parameters
                                         .get_ground_excited_weights(trap_places, &temperature)?;
    

        if transitions.get_delocalised(){
                let mut delocalised_rates: Vec<Candidate> = Vec::with_capacity(trap_places.available_count()*2);
                let inputs = trap_parameters.get_delocalised(
                                                                        &ground_weights, 
                                                                        &excited_weights, 
                                                                        trap_places, 
                                                                        temperature)?;
                delocalised_candidates(&transitions.delocalised, 
                                        inputs, 
                                        trap_places.available(),
                                        &mut delocalised_rates)?;

            
                shortest.find_shortest_from_summed_rates(delocalised_rates, rng)?;
        }
        let (ground_weights, excited_weights) = if trap_places.available_count() > 1 && ground_weights.len() == 1 {
            (vec![ground_weights[0]; trap_places.available_count()], vec![excited_weights[0]; trap_places.available_count()])
        }else{
            (ground_weights, excited_weights)
        };
    
        if transitions.get_localised_recombination(){
            let mut localised_rates: Vec<Candidate> = Vec::with_capacity(trap_places.available_count()*hole_places.available_count()*2);
            let mut distance: Vec<Float> = Vec::with_capacity(hole_places.available_count());
            
            let hole_coords: Vec<Coord> = hole_places.available_indices_vec()
                                                    .into_iter()
                                                    .map(|index|places.holes()[index])
                                                    .collect();
            
            for ((&source, ground_weight),excited_weight) in trap_places.available()
                                                                                        .iter()
                                                                                        .zip(&ground_weights)
                                                                                        .zip(&excited_weights) 
            {
                distance.clear();

                let source_coord = &places.traps()[source.index()];
                distance.extend(hole_coords
                                    .iter()
                                    .map(|hole_coord| 
                                        cube.distance(source_coord, hole_coord)),
                                );
                let inputs = trap_parameters.get_localised_recombination(
                                                                        &ground_weight, 
                                                                        &excited_weight, 
                                                                        &distance, 
                                                                        &source, 
                                                                        temperature)?;
                        
                localised_recombination_candidates(&transitions.localised_recomb, 
                                                inputs, 
                                                source, 
                                                hole_places.available(), 
                                                &mut localised_rates)?;  
            }
    
            shortest.find_shortest_from_summed_rates(localised_rates, rng)?;
        }

        if transitions.get_localised_retrapping(){

            let mut localised_rates: Vec<Candidate> = Vec::with_capacity(trap_places.available_count()*trap_places.unavailable_count()*2);
            
            let mut distance: Vec<Float> = Vec::with_capacity(trap_places.unavailable_count());
            let destination_coords: Vec<Coord> = trap_places.unavailable_indices_vec()
                                                            .into_iter()
                                                            .map(|index|places.traps()[index])
                                                            .collect();

            for ((&source, ground_weight),excited_weight) in trap_places.available()
                                                                                        .iter()
                                                                                        .zip(&ground_weights)
                                                                                        .zip(&excited_weights) {

                distance.clear();            
                let source_coord = &places.traps()[source.index()];
                distance.extend(destination_coords
                                    .iter()
                                        .map(|dcoord| 
                                            cube.distance(source_coord, dcoord)),
                                );
                let inputs = trap_parameters.get_localised_retrapping(
                                                                            &ground_weight, 
                                                                            &excited_weight, 
                                                                            &distance, 
                                                                            &source, 
                                                                            temperature)?;
                
                localised_retrapping_candidates(&transitions.localised_retrap, 
                                                inputs, 
                                                source, 
                                                trap_places.unavailable(), 
                                                &mut localised_rates)?;
            }
            shortest.find_shortest_from_summed_rates(localised_rates, rng)?;
        }
    }

    if transitions.get_filling(){
        filling_inputs.update_occ(trap_places.available_count());
        
        let fill = filling_candidate(
            &transitions.filling,
            &filling_inputs,
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
    retrapping_parameters: &ReTrapParameterLayout,
    cube: &Cube,
    rng: &mut impl Rng,
) -> Result<Event, String> {
    match event {
        Event::LocalisedRecombination { source, hole, .. }
        | Event::DelocalisedRecombination { source, hole, .. } 
        | Event::FillingLoss { trap: source, hole } => {
            if !trap_places.make_unavailable(source) {
                return Err(format!("recombination source {source:?} was not occupied"));
            }
            if !hole_places.make_unavailable(hole) {
                return Err(format!("recombination hole {hole:?} was not available"));
            }
        }
        Event::LocalisedRetrapping { source, destination, .. }
        | Event::DelocalisedRetrapping { source, destination, .. } 
        | Event::FillingTrapOnly { trap_lost: source, trap_gain: destination } => {
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
            let delocal_outcome = choose_delocalised_outcome(
                source,
                places,
                trap_places,
                hole_places,
                retrapping_parameters,
                cube,
                state,
                rng
            )?;
            return Ok(delocal_outcome);   
        }

        Event::FillingSelect { .. } => {
            let filling = choose_filling_outcome(places, trap_places, hole_places, retrapping_parameters, cube, rng)?;
            return Ok(filling)
        }
        Event::FillingStandard{trap, hole} => {
            if !hole_places.make_available(hole) {
                return Err(format!("retrapping hole {hole:?} was not available"));
            }
            if !trap_places.make_available(trap) {
                return Err(format!("filling destination {trap:?} was already occupied"));
            }
        },
        Event::FillingHoleOnly { hole_lost, hole_gain } => {
            if !hole_places.make_unavailable(hole_lost) {
                return Err(format!("recombination source {hole_lost:?} was not occupied"));
            }
            if !hole_places.make_available(hole_gain) {
                return Err(format!("recombination hole {hole_gain:?} was not available"));
            }
        },
    
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
    trap_places: &mut PlaceAvailability,
    hole_places: &mut PlaceAvailability,
    parameters: &ReTrapParameterLayout,
    cube: &Cube,
    state: ElectronicState,
    rng: &mut impl Rng,
)  -> Result<Event , String,> {
        
        if !trap_places.make_unavailable(source) {
                return Err(format!("Delocalised source {source:?} was not occupied"));
        }

        match parameters {
            ReTrapParameterLayout::GaussianReTrappingCBFill { cb_hole_to_trap, cb_mu, ..} | 
            ReTrapParameterLayout::GaussianFillReTrappingCB { cb_hole_to_trap, cb_mu, ..} |
            ReTrapParameterLayout::GaussianReTrappingCB { cb_hole_to_trap, mu: cb_mu }
            => {
                let source_position = &places.traps()[source.index()];
                let mut destinations: Vec<ReciprocalCandidate> = Vec::with_capacity(hole_places.available_count()+trap_places.unavailable_count());
                for &hole in hole_places.available() {
                    let distance = cube.distance(source_position, &places.holes()[hole.index()]);
                    destinations.push(
                        gaussian_kernel(
                            1.0, 
                            *cb_mu, distance, 
                            hole, 
                            true)?);
                }
                for &trap in trap_places.unavailable() {
                let distance = cube.distance(source_position, &places.traps()[trap.index()]);
                destinations.push(
                    gaussian_kernel(
                        *cb_hole_to_trap, 
                        *cb_mu, distance, 
                        trap, 
                        false)?);
                }
                let choice = find_shortest_from_summed_reciprocal_rates(destinations, rng)?;
                if choice.hole {
                    if !hole_places.make_unavailable(choice.place) {
                        return Err(format!(
                            "delocalised recombination hole {:?} was not available", choice.place.index()
                        ));
                    }
                    return Ok(Event::DelocalisedRecombination { source, hole: choice.place, state });

                }else {
                    if !trap_places.make_available(choice.place) {
                        return Err(format!(
                            "delocalised retrapping destination {:?} was occupied", choice.place.index()
                        ));
                    }
                    return Ok(Event::DelocalisedRetrapping { source, destination: choice.place, state  } );

                }
            },
            ReTrapParameterLayout::NoneGaussianReTrappingCBFill { cb_hole_to_trap,.. } |
            ReTrapParameterLayout::NoneGaussianFillReTrappingCB { cb_hole_to_trap } | 
            ReTrapParameterLayout::NoneGaussianReTrappingCB { cb_hole_to_trap }
            => {
                if rng.random_bool(*cb_hole_to_trap){
                    let trap = {
                        let trap_dest = trap_places.unavailable();
                        if trap_dest.is_empty() {
                            return Err("filling selected when no empty traps remain".to_string());
                        }
                        trap_dest[rng.random_range(0..trap_dest.len())]
                    };
                    if !trap_places.make_available(trap) {
                        return Err(format!(
                            "delocalised retrapping destination {:?} was occupied", trap.index()
                        ));
                    }
                    return Ok(Event::DelocalisedRetrapping { source, destination: trap, state });

                } else {
                    let hole = {
                    let empty_holes = hole_places.available();
                    if empty_holes.is_empty() {
                        return Err("No holes to put electron in".to_string());
                    }
                    empty_holes[rng.random_range(0..empty_holes.len())]
                    };
                    if !hole_places.make_unavailable(hole) {
                        return Err(format!(
                            "delocalised recombination hole {:?} was not available", hole.index()
                        ));
                    }
                    return Ok( Event::DelocalisedRecombination { source, hole, state });    
                }
            },
            ReTrapParameterLayout::GaussianCBReTrappingFill {cb_mu, .. } |
            ReTrapParameterLayout::GaussianCBFill { cb_mu, .. } | 
            ReTrapParameterLayout::GaussianCB { mu: cb_mu } 
            => {
          
                let source_position = &places.traps()[source.index()];
                let mut destinations: Vec<ReciprocalCandidate> = Vec::with_capacity(hole_places.available_count()+trap_places.unavailable_count());
                for &hole in hole_places.available() {
                    let distance = cube.distance(source_position, &places.holes()[hole.index()]);
                    destinations.push(
                        gaussian_kernel(
                            1.0, 
                            *cb_mu, distance, 
                            hole, 
                            true)?);
                }
                let choice = find_shortest_from_summed_reciprocal_rates(destinations, rng)?;
                if !hole_places.make_unavailable(choice.place) {
                        return Err(format!(
                            "delocalised recombination hole {:?} was not available", choice.place.index()
                        ));
                    }
                return Ok(Event::DelocalisedRecombination { source, hole: choice.place, state  });

            },
            ReTrapParameterLayout::NoneGaussianCBReTrappingFill { .. } |
            ReTrapParameterLayout::NoneGaussianCBFill |
            ReTrapParameterLayout::NoneGaussianCB => {
                let hole = {
                let empty_holes = hole_places.available();
                if empty_holes.is_empty() {
                    return Err("No holes to put electron in".to_string());
                }
                empty_holes[rng.random_range(0..empty_holes.len())]
                };
                if !hole_places.make_unavailable(hole) {
                        return Err(format!(
                            "delocalised recombination hole {:?} was not available", hole.index()
                        ));
                    }
                return Ok( Event::DelocalisedRecombination { source, hole, state },);
            },
            _ => {
                return Err(format!("A delocalsied transition destination has been asked for but delocalised transitions are not on. 
                                    This should never happen!"))
            }

        }

    }


/// Resolve an aggregate irradiation event into population changes.
///
/// Irradiation first activates one previously inactive hole site. Normally it
/// also occupies a uniformly selected empty trap. When filling-time
/// recombination is enabled, a fixed 0.5 branch instead consumes an active
/// hole; the returned [`Event::Filling`] records the selected identifiers.
pub fn choose_filling_outcome(
    places: &ElectronPlaces,
    trap_places: &mut PlaceAvailability,
    hole_places: &mut PlaceAvailability,
    parameters: &ReTrapParameterLayout,
    cube: &Cube,
    rng: &mut impl Rng,
) -> Result<Event, String> {
    
    match parameters {
        ReTrapParameterLayout::GaussianCBReTrappingFill{cb_hole_to_trap,cb_mu,vb_trap_to_hole, vb_mu} |
        ReTrapParameterLayout::GaussianReTrappingFill{cb_hole_to_trap,cb_mu,vb_trap_to_hole, vb_mu} |
        ReTrapParameterLayout::GaussianReTrappingCBFill{cb_hole_to_trap,cb_mu,vb_trap_to_hole, vb_mu} 
        => {
            let source_position = Coord::random_in(cube.boundary.x, cube.boundary.y, cube.boundary.z, rng)?;
            
            let mut destinations: Vec<ReciprocalCandidate> = Vec::with_capacity(hole_places.unavailable_count()+trap_places.available_count());
            for &hole in hole_places.unavailable() {
                let distance = cube.distance(&source_position, &places.holes()[hole.index()]);
                destinations.push(
                    gaussian_kernel(
                        1.0, 
                        *vb_mu, distance, 
                        hole,
                    true)?);
            }
            for &trap in trap_places.available() {
                let distance = cube.distance(&source_position, &places.traps()[trap.index()]);
                destinations.push(
                    gaussian_kernel(
                        *vb_trap_to_hole, 
                        *vb_mu, distance, 
                        trap, false)?);
            }

            let hole = find_shortest_from_summed_reciprocal_rates(destinations, rng)?;
            if hole.hole {
                if !hole_places.make_available(hole.place) {
                    return Err(format!("filling destination {:?} was occupied",hole.place));
                }  
            }else {
                if !trap_places.make_unavailable(hole.place) {
                    return Err(format!("filling destination {:?} was occupied",hole.place));
                }
            }
            
            let mut destinations: Vec<ReciprocalCandidate> = Vec::with_capacity(trap_places.unavailable_count()+hole_places.available_count());
            for &trap in trap_places.unavailable() {
                let distance = cube.distance(&source_position, &places.traps()[trap.index()]);
                destinations.push(
                    gaussian_kernel(
                        1.0, 
                        *cb_mu, distance, 
                        trap, false)?);
            }
            for &hole in hole_places.available() {
                let distance = cube.distance(&source_position, &places.holes()[hole.index()]);
                destinations.push(
                    gaussian_kernel(
                        *cb_hole_to_trap, 
                        *cb_mu, distance, 
                        hole,
                    true)?);
            }
            let trap = find_shortest_from_summed_reciprocal_rates(destinations, rng)?;
            if trap.hole {
                if !hole_places.make_unavailable(trap.place) {
                    return Err(format!("filling destination {:?} was occupied",trap.place));
                }
            } else {
                if !trap_places.make_available(trap.place) {
                    return Err(format!("filling destination {:?} was occupied",trap.place));
                }
                
            }
            match (hole.hole, trap.hole){
                (true, true)   => return Ok(Event::FillingHoleOnly { hole_lost: trap.place, hole_gain: hole.place }),
                (true, false)  => return Ok(Event::FillingStandard { trap: trap.place, hole: hole.place }),
                (false, true)  => return Ok(Event::FillingLoss { trap: hole.place, hole:trap.place }),
                (false, false) => return Ok(Event::FillingTrapOnly { trap_lost: hole.place, trap_gain: trap.place }),
            }
            
        },
        ReTrapParameterLayout::NoneGaussianCBReTrappingFill{cb_hole_to_trap,vb_trap_to_hole} |
        ReTrapParameterLayout::NoneGaussianReTrappingFill{cb_hole_to_trap,vb_trap_to_hole} |
        ReTrapParameterLayout::NoneGaussianReTrappingCBFill{cb_hole_to_trap,vb_trap_to_hole}
        => {
            // If true electron in conduction band goes straight to a hole
            let cb_hole = rng.random_bool(*cb_hole_to_trap);
            // If true the hole has gone to an occupied trap and annihilated it
            let vb_trap = rng.random_bool(*vb_trap_to_hole);

            match (cb_hole, vb_trap){
            (true, true) => { // electron to hole, hole to trap
                let hole = {
                    let empty_holes = trap_places.available();
                    if empty_holes.is_empty() {
                        return Ok(Event::None);
                        // return Err("Filling selected but no traps for hole to annihilate remain".to_string());
                    }
                    empty_holes[rng.random_range(0..empty_holes.len())]
                };

                if !trap_places.make_unavailable(hole) {
                    return Err(format!("Hole annihilation trap {hole:?} was not occupied"));
                }

                let trap = {
                    let trap_dest = hole_places.available();
                    if trap_dest.is_empty() {
                        trap_places.make_available(hole);
                        return Ok(Event::None);
                        // return Err("Filling selected but no holes for electron to recombine with remain".to_string());
                    }
                    trap_dest[rng.random_range(0..trap_dest.len())]
                };

                if !hole_places.make_unavailable(trap) {
                    return Err(format!("Hole to be recombined with via filling {trap:?} was occupied"));
                }
                return Ok(Event::FillingLoss { trap, hole });
            },
            (true, false) => { // electron to hole, hole to hole
                let hole = {
                    let empty_holes = hole_places.unavailable();
                    if empty_holes.is_empty() {
                        return Ok(Event::None);
                        // return Err("Filling selected when no available holes remain".to_string());
                    }
                    empty_holes[rng.random_range(0..empty_holes.len())]
                };

                if !hole_places.make_available(hole) {
                    return Err(format!("Filling destination {hole:?} was occupied"));
                }
                let trap = {
                    let trap_dest = hole_places.available();
                    if trap_dest.is_empty() {
                        hole_places.make_unavailable(hole);
                        return Ok(Event::None);
                        // return Err("Filling recombination selected when no available holes remaining".to_string());
                    }
                    trap_dest[rng.random_range(0..trap_dest.len())]
                };

                if !hole_places.make_unavailable(trap) {
                    return Err(format!("Filling recombination destination {trap:?} was already occupied"));
                }
                return Ok(Event::FillingHoleOnly { hole_lost: trap, hole_gain: hole });
            },
            (false, true) => { // electron to trap, hole to trap
                let trap = {
                    let trap_dest = trap_places.unavailable();
                    if trap_dest.is_empty() {
                        return Ok(Event::None);
                        // return Err("Filling selected when no empty traps remain".to_string());
                    }
                    trap_dest[rng.random_range(0..trap_dest.len())]
                };

                if !trap_places.make_available(trap) {
                    return Err(format!("Filling destination {trap:?} was occupied"));
                }

                let hole = {
                    let empty_holes = trap_places.available();
                    if empty_holes.is_empty() {
                        trap_places.make_unavailable(trap);
                        return Ok(Event::None);

                        // return Err("Filling trap annihilation selected when no traps are occupied".to_string());
                    }
                    empty_holes[rng.random_range(0..empty_holes.len())]
                };

                if !trap_places.make_unavailable(hole) {
                    return Err(format!("Filling trap annihilation destination {hole:?} was empty"));
                }
                
                return Ok(Event::FillingTrapOnly { trap_lost: hole, trap_gain: trap });
            },
            (false, false) =>  { // electron to trap, hole to hole
                let hole = {
                    let empty_holes = hole_places.unavailable();
                    if empty_holes.is_empty() {
                        return Ok(Event::None);
                        // return Err("Filling selected when no available holes remain".to_string());
                    }
                    empty_holes[rng.random_range(0..empty_holes.len())]
                };

                if !hole_places.make_available(hole) {
                    return Err(format!("Filling destination {hole:?} was occupied"));
                }
                let trap = {
                    let trap_dest = trap_places.unavailable();
                    if trap_dest.is_empty() {
                        hole_places.make_unavailable(hole);
                        return Ok(Event::None);
                        // return Err("Filling selected when no empty traps remain".to_string());
                    }
                    trap_dest[rng.random_range(0..trap_dest.len())]
                };

                if !trap_places.make_available(trap) {
                    return Err(format!("Filling destination {trap:?} was occupied"));
                }
                return Ok(Event::FillingStandard { trap, hole });
            }
        }

        },
        ReTrapParameterLayout::GaussianFill{cb_mu,vb_mu} |
        ReTrapParameterLayout::GaussianCBFill{cb_mu, vb_mu} |
        ReTrapParameterLayout::GaussianFillReTrappingCB{cb_mu, vb_mu, ..} 
        => {
            let source_position = Coord::random_in(cube.boundary.x, cube.boundary.y, cube.boundary.z, rng)?;
            let mut destinations: Vec<ReciprocalCandidate> = Vec::with_capacity(hole_places.unavailable_count());
                for &hole in hole_places.unavailable() {
                    let distance = cube.distance(&source_position, &places.holes()[hole.index()]);
                    destinations.push(
                        gaussian_kernel(
                            1.0, 
                            *vb_mu, distance, 
                            hole,
                        true)?);
                }
            let hole = find_shortest_from_summed_reciprocal_rates(destinations, rng)?;
            if !hole_places.make_available(hole.place) {
                return Err(format!("Filling destination {:?} was occupied",hole.place.index()));
            }
            
            let mut destinations: Vec<ReciprocalCandidate> = Vec::with_capacity(trap_places.unavailable_count());
                for &trap in trap_places.unavailable() {
                    let distance = cube.distance(&source_position, &places.traps()[trap.index()]);
                    destinations.push(
                        gaussian_kernel(
                            1.0, 
                            *cb_mu, distance, 
                            trap, false)?);
            }
            let trap = find_shortest_from_summed_reciprocal_rates(destinations, rng)?;
            if !trap_places.make_available(trap.place) {
                return Err(format!("Filling destination {:?} was occupied",trap.place.index()));
            }
            return Ok(Event::FillingStandard { trap: trap.place, hole: hole.place,}); 
           
        },

        ReTrapParameterLayout::NoneGaussianFill |
        ReTrapParameterLayout::NoneGaussianCBFill |
        ReTrapParameterLayout::NoneGaussianFillReTrappingCB{..} 
        => {
            let hole = {
            let empty_holes = hole_places.unavailable();
            if empty_holes.is_empty() {
                return Err("Filling selected when no available holes remain".to_string());
            }

            empty_holes[rng.random_range(0..empty_holes.len())]
            };

            if !hole_places.make_available(hole) {
                return Err(format!("Filling destination {hole:?} was occupied"));
            }
            let trap = {
                let trap_dest = trap_places.unavailable();
                if trap_dest.is_empty() {
                    return Err("Filling selected when no empty traps remain".to_string());
                }
                trap_dest[rng.random_range(0..trap_dest.len())]
            };

            if !trap_places.make_available(trap) {
                return Err(format!("Filling destination {trap:?} was occupied"));
            }
            return Ok(Event::FillingStandard { trap, hole, });
        },

        _ => {return Err(format!("A filling event has been asked for but filling transitions are not on. This should never happen! The type was {:?}",parameters))}

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
    retrapping_parameters: &ReTrapParameterLayout,
    filling_inputs: &mut FillingTransitionInputs,
    time_temperature: &mut TimeTemperature,
    cube: &Cube,
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
            filling_inputs,
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
                    retrapping_parameters,
                    cube,
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
