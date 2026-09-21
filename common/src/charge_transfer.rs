// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Physical events and stochastic waiting times used by the kinetic Monte Carlo engine.
//!
//! A [`crate::charge_transfer::Candidate`] couples a possible charge-transfer
//! event to its rate. For a process with constant rate `k`,
//! [`crate::charge_transfer::TimedCandidate::rate_to_lifetime`]
//! samples the exponentially distributed waiting time `tau = -ln(U) / k`.
//! Competing pathways are handled by retaining the candidate with the shortest
//! sampled lifetime.

use crate::place_ids::{PlaceId};
use crate::numeric::{TimeFloat,Float};
use serde::{Deserialize, Serialize};
use crate::trap_hole_band_tail::TrapParameters;
use crate::rate_equation_inputs::{
    DelocalisedTransitionInputs, FillingTransitionInputs, LocalisedTransitionInputs,
};
use crate::rate_equation_selection::{
    DelocalisedRateEquation, FillingRateEquation, LocalisedRateEquation,
};
use crate::random::generatre_exponential_random;
use crate::rate_equations::{retrapping_probability_by_r};
use rand::{Rng, RngExt};



/// Sentinel rate used internally for a transition pathway disabled by configuration.
///
/// Physical rates must be non-negative. [`TimedCandidate::find_shortest`]
/// recognises this value and omits the candidate from the competition.
pub const DISABLED_RATE: TimeFloat = -1.0;

/// Localised state occupied immediately before a charge-transfer event.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElectronicState {
    /// Lower-energy trap state.
    Ground,
    /// Thermally populated excited trap state.
    Excited,
}

/// Discrete state change produced by one kinetic Monte Carlo event.
///
/// Site identifiers refer to the trap and hole collections belonging to the
/// same spatial realisation. Delocalised release is represented separately
/// until its recombination or retrapping destination has been sampled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Event {
    /// An electron tunnels directly from a trap to a hole centre.
    LocalisedRecombination {
        /// Initially occupied electron trap.
        source: PlaceId,
        /// Hole centre consumed by recombination.
        hole: PlaceId,
        /// Trap state from which the electron tunnelled.
        state: ElectronicState,
    },
    /// An electron tunnels directly between two localised traps.
    LocalisedRetrapping {
        /// Initially occupied electron trap.
        source: PlaceId,
        /// Initially empty trap that receives the electron.
        destination: PlaceId,
        /// Trap state from which the electron tunnelled.
        state: ElectronicState,
    },
    /// An electron is thermally released from a trap into the conduction band.
    Delocalised {
        /// Trap that released the electron.
        source: PlaceId,
        /// Trap state from which the electron was released.
        state: ElectronicState,
    },
    /// A conduction-band electron recombines at a hole centre.
    DelocalisedRecombination {
        /// Trap from which the electron entered the conduction band.
        source: PlaceId,
        /// Hole centre consumed by recombination.
        hole: PlaceId,
        /// Original ground or excited state of the released electron.
        state: ElectronicState,
    },
    /// A conduction-band electron is captured by an empty trap.
    DelocalisedRetrapping {
        /// Trap from which the electron entered the conduction band.
        source: PlaceId,
        /// Empty trap that captures the electron.
        destination: PlaceId,
        /// Original ground or excited state of the released electron.
        state: ElectronicState,
    },
    /// Irradiation creates a hole and either fills a trap or recombines a carrier.
    Filling {
        /// Destination identifier stored by the current filling branch.
        ///
        /// This is a trap ID for ordinary filling. The filling-time
        /// recombination branch currently stores its consumed hole ID here.
        trap: PlaceId,
        /// Newly activated hole site selected by the filling model.
        hole: PlaceId,
    },
    /// Profile boundary or initial sample with no charge-transfer event.
    None,
}


/// One possible state transition and its Poisson rate in inverse seconds.
#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    /// State change that occurs if this candidate wins the lifetime competition.
    pub event: Event,
    /// Hazard rate used to sample the candidate lifetime, in s⁻¹.
    pub rate: TimeFloat,
}


/// Build ground- and excited-state thermal-release candidates for one trap.
/// 
/// Delocalised rates are calculated using [crate::rate_equation_selection::DelocalisedRateEquation::calculate]
/// Each Arrhenius release rate is weighted by the thermal occupation of
/// its originating state. A disabled state receives [`DISABLED_RATE`].
pub fn delocalised_candidates<'a>(
    transitions: &DelocalisedRateEquation,
    inputs: DelocalisedTransitionInputs<Vec<Float>, Vec<Float>, &'a [Float]>,
    sources: &[PlaceId],
    output: &mut Vec<Candidate>,
    
) -> Result<(), String>{
    
    let (ground_rate, excited_rate) = transitions.calculate(&inputs);

    let ground_rates = ground_rate
                            .ok_or_else(|| {format!(
        "could not calculate ground-state rate: \
         e_cb={}, frequency={}, weights={}",
        inputs.e_cb_ground.len(),
        inputs.frequency_ground.len(),
        inputs.ground_weight.len(),
    )})?;

    let excited_rates = excited_rate
                            .ok_or_else(|| "could not calculate excited-state rate".to_string())?;
    
    if sources.len() > 1 && ground_rates.len() == 1 {
        for &source in sources {
            output.push(Candidate {
            event: Event::Delocalised {
            source,
            state: ElectronicState::Ground,
        },
        rate: ground_rates[0],
        });

        output.push(Candidate {
            event: Event::Delocalised {
                source,
                state: ElectronicState::Excited,
            },
            rate: excited_rates[0],
        });
        }
    }else{
        for ((&source, ground_rate), excited_rate) in sources
            .iter()
            .zip(ground_rates)
            .zip(excited_rates)
        {  output.push(Candidate {
                event: Event::Delocalised {
                source,
                state: ElectronicState::Ground,
            },
            rate: ground_rate,
            });

            output.push(Candidate {
                event: Event::Delocalised {
                    source,
                    state: ElectronicState::Excited,
                },
                rate: excited_rate,
            });
        }
    }
    Ok(())

}

/// Build ground- and excited-state tunnelling candidates to one hole.
///
/// Localised recombination rates are calculated using [crate::rate_equation_selection::LocalisedRateEquation::calculate]
/// The selected equation applies exponential distance attenuation and the
/// state occupation weights to the trap-to-hole separation.
pub fn localised_recombination_candidates<'a>(
    transitions: &LocalisedRateEquation,
    inputs: LocalisedTransitionInputs<Float, Float, &'a Float, &'a [Float]>,
    source: PlaceId,
    holes: &[PlaceId],
    output: &mut Vec<Candidate>,
) -> Result<(), String> {

    
    let (ground_rate, excited_rate) = transitions.calculate(&inputs);

    let ground_rates = ground_rate
                            .ok_or_else(|| "could not calculate ground-state rate".to_string())?;
    
    let excited_rates = excited_rate
                            .ok_or_else(|| "could not calculate excited-state rate".to_string())?;
    
    for ((&hole, ground_rate), excited_rate) in holes
        .iter()
        .zip(ground_rates)
        .zip(excited_rates)
    { 
        output.push(Candidate {
            event: Event::LocalisedRecombination {
            source,
            hole,
            state: ElectronicState::Ground,
        },
        rate: ground_rate,
        });

        output.push(Candidate {
            event: Event::LocalisedRecombination {
                source,
                hole,
                state: ElectronicState::Excited,
            },
            rate: excited_rate,
        });
    }
    
    Ok(())

}

/// Build ground- and excited-state tunnelling candidates to an empty trap.
///
/// Localised retrapping rates are calculated using [crate::rate_equation_selection::LocalisedRateEquation::calculate]
/// The rate has the same distance-decay form as localised recombination,
/// but its event transfers occupancy between trap identifiers.
pub fn localised_retrapping_candidates<'a>(
    transitions: &LocalisedRateEquation,
    inputs: LocalisedTransitionInputs<Float, Float, &'a Float, &'a [Float]>,
    source: PlaceId,
    destinations: &[PlaceId],
    output: &mut Vec<Candidate>,
) -> Result<(), String>{

    let (ground_rate, excited_rate) = transitions.calculate(&inputs);
    
    let ground_rates = ground_rate
                            .ok_or_else(|| "could not calculate ground-state rate".to_string())?;
    let excited_rates = excited_rate
                            .ok_or_else(|| "could not calculate excited-state rate".to_string())?;

    for ((&destination, ground_rate), excited_rate) in destinations
        .iter()
        .zip(ground_rates)
        .zip(excited_rates)
    { 
        output.push(Candidate {
            event: Event::LocalisedRetrapping {
            source,
            destination,
            state: ElectronicState::Ground,
        },
        rate: ground_rate,
        });

        output.push(Candidate {
            event: Event::LocalisedRetrapping {
                source,
                destination,
                state: ElectronicState::Excited,
            },
            rate: excited_rate,
        });
    }
    
    Ok(())
    
}

impl Candidate{
    /// Build the aggregate dose-driven filling candidate for all empty traps.
    ///
    /// Filling rates are calculated using [crate::rate_equation_selection::FillingRateEquation::calculate]
    /// The underlying rate is `(dose_rate / characteristic_dose) *
    /// (total_population - occupied_population)`. Concrete electron and hole
    /// destinations are chosen only if this aggregate candidate fires.
    pub fn filling_candidate(
        transitions: &FillingRateEquation,
        d0: &Float, d_dot: &Float,
        occupied_population: usize,
        total_population: usize,  
    ) -> Result<Self, String> {

        let filling_inputs = FillingTransitionInputs {
            characteristic_dose: *d0,
            dose_rate: *d_dot,
            occupied_population: occupied_population as Float,
            total_population: total_population as Float,
        };

        let rate: Option<TimeFloat> = transitions.calculate(&filling_inputs);
        
        let rate = rate
                        .ok_or_else(|| "could not calculate filling rate".to_string())?;
   
        Ok(
            Self{ event: Event::Filling { trap: PlaceId::new(0)?, hole: PlaceId::new(0)?}, rate: rate}
        )
    
    }
        /// Sample the destination time for conduction-band recombination.
    ///
    /// The distance model supplies a reciprocal rate proportional to
    /// `prefactor * exp((distance / mu)^2)`. Multiplying this by `-ln(U)`
    /// therefore makes nearby centres statistically more likely to win.
    pub fn delocalised_recombination(
        prefactor: Float,
        mu: Float,
        distance: Float,
        source: PlaceId,
        hole: PlaceId,
        state: ElectronicState,
    ) -> Result<Self, String> {
        
        let reciprocal_rate: Option<TimeFloat> = retrapping_probability_by_r(&prefactor, &mu, &distance);
        let reciprocal_rate = reciprocal_rate.ok_or_else(|| {
            format!("could not calculate delocalised destination rate for hole: {}", hole.index())
        })?;

        if reciprocal_rate.is_nan() || reciprocal_rate < 0.0 || reciprocal_rate.is_infinite(){
             return Ok(Self  { event: Event::DelocalisedRecombination { source, hole, state  }, rate: DISABLED_RATE});
        }
       
        return Ok(Candidate  { event: Event::DelocalisedRecombination { source, hole, state  }, rate: 1.0/reciprocal_rate});
        
    }


}

/// A candidate after its stochastic waiting time has been sampled.
#[derive(Debug, Clone, Copy)]
pub struct TimedCandidate {
    /// State change associated with the sampled lifetime.
    pub event: Event,
    /// Non-negative waiting time from the current simulation state, in seconds.
    pub time: TimeFloat,
}

impl TimedCandidate {
    /// Create an inert candidate at infinite time.
    ///
    /// Despite the historical method name, the returned time is positive
    /// infinity. It is useful as the initial value of a minimum search.
    pub fn new_negative_time() -> Self {
        Self { event: Event::None, time: TimeFloat::INFINITY,}
    }

    /// Sample an exponentially distributed lifetime for a constant-rate event.
    ///
    /// The inverse-transform expression is `tau = -ln(U) / rate`, with `U`
    /// uniformly distributed on the open interval `(0, 1)`. A zero rate maps
    /// to infinite time; negative and non-finite rates are rejected.
    pub fn rate_to_lifetime(candidate : Candidate, rng: &mut impl Rng) -> Result<Self, String> {
        if !candidate.rate.is_finite() {
            return Err(format!("non-finite transition rate: {}", candidate.rate));
        }

        if candidate.rate < 0.0 {
            return Err(format!("negative transition rate: {}", candidate.rate));
        }

        if candidate.rate == 0.0 {
            return Ok(TimedCandidate { event: candidate.event, time: TimeFloat::INFINITY });
        }
        let u = generatre_exponential_random(rng)as TimeFloat;
        
        // let u: TimeFloat = rng.sample(rand::distributions::Open01);
        let time = u / candidate.rate;
        Ok(TimedCandidate { event: candidate.event, time})

    }

    /// Sample `candidate` and retain it when its lifetime is shorter than `self`.
    ///
    /// A rate equal to [`DISABLED_RATE`] is skipped without consuming random
    /// numbers, keeping disabled pathways outside the stochastic competition.
    pub fn find_shortest(&mut self, candidate: Candidate, rng: &mut impl Rng) -> Result<(), String> {
        if candidate.rate == DISABLED_RATE {
            return Ok(());
        }
        let tc = TimedCandidate::rate_to_lifetime(candidate,rng)?;
        if tc.time.total_cmp(&self.time).is_lt(){
            *self = tc;
            return Ok(());
        }else {
            return Ok(());
        }

    }
    /// Takes a vector of candidate weights are summed, a rate equal to [`DISABLED_RATE`], or is not finite or negative is skipped.
    ///                                             r0 =∑r_i
    ///                                         τ = -ln(u1)/r0
    /// A trial time is then found with a single random number call. If this time is less than the current smallest time
    /// a second random number is generated and multiplied by the total rate. A specific event is chosen according to 
    ///                 ∑^p_i=1 r_i < u2*r0 ≤ ∑^p _i=1 r_i

    pub fn find_shortest_from_summed_rates(&mut self, candidates: Vec<Candidate>, rng: &mut impl Rng) -> Result<(), String> {
        
        let total_rate: TimeFloat = candidates
        .iter()
        .filter(|candidate| {
                candidate.rate != DISABLED_RATE
                && candidate.rate.is_finite()
                && candidate.rate >= 0.0
            })
        .map(|candidate| candidate.rate)
        .sum();

        if total_rate == 0.0 {
            return Ok(());
        }
        let u = generatre_exponential_random(rng)as TimeFloat;
        let event_time: TimeFloat = u / total_rate ;
     
        if event_time.total_cmp(&self.time).is_lt(){
            
            let u_event: TimeFloat = rng.sample(rand::distr::Open01);
            let mut target = (u_event * total_rate) as TimeFloat;
            
            for candidate in candidates {
                if candidate.rate <= 0.0 {
                    continue;
                }

                if target < candidate.rate {
                    *self = TimedCandidate { event: candidate.event, time: event_time};
                    return Ok(());
                }

                target -= candidate.rate;
            }

            return Ok(());
            
          
        }else {
            return Ok(());
        }

    }

    /// Retain an already sampled candidate when it occurs sooner than `self`.
    pub fn find_smallest_candidate(&mut self, candidate: TimedCandidate)-> Result<(), String> { 
        if candidate.time.total_cmp(&self.time).is_lt(){
            *self = candidate;
            return Ok(());
        }else {
            return Ok(());
        }

    }
    
}

/// Serializable snapshot emitted after a profile step or physical event.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RecordedEvent {
    /// Absolute profile time in seconds.
    pub time: TimeFloat,
    /// Fraction of traps occupied by electrons after the event.
    pub fill: Float,
    /// Profile temperature in kelvin after the event.
    pub temperature: Float,
    /// Physical transition applied at this time, or [`Event::None`] at a boundary.
    pub event: Event,
}
