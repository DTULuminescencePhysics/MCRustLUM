// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use crate::place_ids::{PlaceId};
use crate::numeric::{TimeFloat,Float};
use serde::{Deserialize, Serialize};
use rand::Rng;

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
pub enum DelocalisedOutcome {
    Recombination { hole: PlaceId },
    Retrapping { destination: PlaceId },
}

#[derive(Debug, Clone, Copy)]
pub struct TimedDelocalisedOutcome {
    pub outcome: DelocalisedOutcome,
    pub time: TimeFloat,
}

#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    pub event: Event,
    pub rate: TimeFloat,
}

#[derive(Debug, Clone, Copy)]
pub struct TimedCandidate {
    pub event: Event,
    pub time: TimeFloat,
}

impl Candidate{

    pub fn lifetime<R: Rng + ?Sized>(rate: TimeFloat, rng: &mut R) -> Result<TimeFloat, String> {
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


}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RecordedEvent {
    pub time: TimeFloat,
    pub fill: Float,
    pub temperature: Float,
    pub event: Event,
}