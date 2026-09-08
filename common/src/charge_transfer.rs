// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

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
use crate::rate_equations::{retrapping_probability_by_r};
use rand::Rng;


pub const DISABLED_RATE: TimeFloat = -1.0;

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


#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    pub event: Event,
    pub rate: TimeFloat,
}

impl Candidate{
   
    pub fn delocalised_candidates(
        transitions: &DelocalisedRateEquation,
        parameters: &TrapParameters,
        source: PlaceId,
        temperature: Float,
        ground_weight: Float,
        excited_weight: Float,
    ) -> Result<(Self, Self), String>{

        let inputs = DelocalisedTransitionInputs {
            e_cb_ground: parameters.e_cb_ground,
            frequency_ground: parameters.de_frequency_ground,
            e_cb_excited: parameters.e_cb_excited,
            frequency_excited: parameters.de_frequency_excited,
            temperature,
            ground_weight,
            excited_weight,
        };
        
        let (ground_rate, excited_rate) = transitions.calculate(&inputs);

        let ground_rate = ground_rate
                               .ok_or_else(|| "could not calculate ground-state rate".to_string())?;

        let excited_rate = excited_rate
                                .ok_or_else(|| "could not calculate excited-state rate".to_string())?;
       
        Ok((
            Self {event: Event::Delocalised { source, state: ElectronicState::Ground}, rate: ground_rate},
            Self {event: Event::Delocalised { source, state: ElectronicState::Excited}, rate: excited_rate}
        ))
    
    }

    pub fn localised_recombination_candidates(
        transitions: &LocalisedRateEquation,
        parameters: &TrapParameters,
        source: PlaceId,
        hole: PlaceId,
        temperature: Float,
        distance: Float,
        ground_weight: Float,
        excited_weight: Float,
    ) -> Result<(Self, Self), String>{

        let inputs = LocalisedTransitionInputs {
            alpha_ground: parameters.alpha_ground,
            frequency_ground: parameters.lo_frequency_ground,
            alpha_excited: parameters.alpha_excited,
            frequency_excited: parameters.lo_frequency_excited,
            ground_weight,
            excited_weight,
            distance,
        };

        let (ground_rate, excited_rate) = transitions.calculate(&inputs);

        let ground_rate = ground_rate
                               .ok_or_else(|| "could not calculate ground-state rate".to_string())?;
        
        let excited_rate = excited_rate
                                .ok_or_else(|| "could not calculate excited-state rate".to_string())?;
        
        Ok((
            Self {event: Event::LocalisedRecombination { 
                                source, 
                                hole, 
                                state: ElectronicState::Ground }, 
                 rate: ground_rate},
            Self {event: Event::LocalisedRecombination { 
                                source, 
                                hole, 
                                state: ElectronicState::Excited },
                  rate: excited_rate}
        ))
    
    }

    pub fn localised_retrapping_candidates(
        transitions: &LocalisedRateEquation,
        parameters: &TrapParameters,
        source: PlaceId,
        destination: PlaceId,
        temperature: Float,
        distance: Float,
        ground_weight: Float,
        excited_weight: Float,
    ) -> Result<(Self, Self), String>{

        let inputs = LocalisedTransitionInputs {
            alpha_ground: parameters.alpha_ground,
            frequency_ground: parameters.lo_frequency_ground,
            alpha_excited: parameters.alpha_excited,
            frequency_excited: parameters.lo_frequency_excited,
            ground_weight,
            excited_weight,
            distance,
        };

        let (ground_rate, excited_rate) = transitions.calculate(&inputs);
       
        let ground_rate = ground_rate
                               .ok_or_else(|| "could not calculate ground-state rate".to_string())?;
        let excited_rate = excited_rate
                                .ok_or_else(|| "could not calculate excited-state rate".to_string())?;

        Ok((
            Self {event: Event::LocalisedRetrapping { 
                                source, 
                                destination, 
                                state: ElectronicState::Ground
                            }, 
                 rate: ground_rate},
            Self { event: Event::LocalisedRetrapping { 
                                 source, 
                                 destination, 
                                 state: ElectronicState::Excited },
                  rate: excited_rate}
        ))
    
    }

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

}
#[derive(Debug, Clone, Copy)]
pub struct TimedCandidate {
    pub event: Event,
    pub time: TimeFloat,
}

impl TimedCandidate {

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

        let u: TimeFloat = rng.sample(rand::distributions::Open01);
        let time = -u.ln() / candidate.rate;
        Ok(TimedCandidate { event: candidate.event, time})

    }

    pub fn delocalised_recombination(
        prefactor: Float,
        mu: Float,
        distance: Float,
        source: PlaceId,
        hole: PlaceId,
        state: ElectronicState,
        rng: &mut impl Rng,
    ) -> Result<Self, String> {
        
        let reciprocal_rate: Option<TimeFloat> = retrapping_probability_by_r(&prefactor, &mu, &distance);
        let reciprocal_rate = reciprocal_rate.ok_or_else(|| {
            format!("could not calculate delocalised destination rate for hole: {}", hole.index())
        })?;

        if reciprocal_rate.is_nan() || reciprocal_rate < 0.0 {
            return Err(format!(
                "invalid delocalised destination reciprocal rate: {reciprocal_rate}"
            ));
        }
        if reciprocal_rate.is_infinite() {
            return Ok(TimedCandidate  { event: Event::DelocalisedRecombination { source, hole, state  }, time: TimeFloat::INFINITY});
        }

        let u: TimeFloat = rng.sample(rand::distributions::Open01);
        let time = -u.ln() * reciprocal_rate;

        return Ok(TimedCandidate  { event: Event::DelocalisedRecombination { source, hole, state  }, time});
        
    }
    pub fn delocalised_retrapping(
        prefactor: Float,
        mu: Float,
        distance: Float,
        source: PlaceId,
        destination: PlaceId,
        state: ElectronicState,
        rng: &mut impl Rng,
    ) -> Result<Self, String> {
        
        let reciprocal_rate: Option<TimeFloat> = retrapping_probability_by_r(&prefactor, &mu, &distance);
        let reciprocal_rate = reciprocal_rate.ok_or_else(|| {
            format!("could not calculate delocalised destination rate for trap: {}", destination.index())
        })?;

        if reciprocal_rate.is_nan() || reciprocal_rate < 0.0 {
            return Err(format!(
                "invalid delocalised destination reciprocal rate: {reciprocal_rate}"
            ));
        }
        if reciprocal_rate.is_infinite() {
            return Ok(TimedCandidate  { event: Event::DelocalisedRetrapping { source, destination, state }, time: TimeFloat::INFINITY});
        }

        let u: TimeFloat = rng.sample(rand::distributions::Open01);
        let time = -u.ln() * reciprocal_rate;

        return Ok(TimedCandidate  { event: Event::DelocalisedRetrapping { source, destination, state  }, time});
        
    }



}

pub struct TimedCandidatePool {
    candidates: Vec<TimedCandidate>,
}
impl TimedCandidatePool {

    pub fn new() -> Self{ 
        Self {
            candidates: Vec::new(),
        }
    }  

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            candidates: Vec::with_capacity(capacity),
        }
    }
    pub fn push_to_lifetime(&mut self, candidate : Candidate, rng: &mut impl Rng) -> Result<(), String> {
        if candidate.rate == DISABLED_RATE {
            return Ok(());
        }
        self.push(TimedCandidate::rate_to_lifetime(candidate,rng)?);
        Ok(())
       
    }
    pub fn push(&mut self, to_push: TimedCandidate){
        self.candidates.push(to_push)
    }

    pub fn earliest_candidate(&self,) -> Option<TimedCandidate> {
        self.candidates
        .iter()
        .copied()
        .min_by(|a, b| a.time.total_cmp(&b.time))
    }


}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RecordedEvent {
    pub time: TimeFloat,
    pub fill: Float,
    pub temperature: Float,
    pub event: Event,
}