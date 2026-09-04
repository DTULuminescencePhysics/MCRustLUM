// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Place Ids are used to uniquely identify each trap, hole or bandtail state 
//! 

use crate::numeric::Float;
use serde::{Deserialize, Serialize};
use rand::Rng;

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
    pub fn set_initial_condition(count: usize, available_count: usize, rng: &mut impl Rng) -> Result<Self, String> {
        let mut places = PlaceAvailability::new(count)?;
        if available_count == 0 {
            return Ok(places);
        } else if available_count == count {
            places.mark_all_available();
            return Ok(places);
        } else {
            places.randomly_make_available(available_count, rng)?;
            return Ok(places);
        }
    }
    pub fn randomly_make_available(&mut self, n: usize, rng: &mut impl Rng) -> Result<(), String> {
        let unavailable_count = self.ids.len() - self.available_count;
        if n > unavailable_count {
            return Err(format!(
                "cannot make {n} places available: only \
                {unavailable_count} places remain"
            ));
        }

        let first_new = self.available_count;
        let new_available_end = first_new + n;
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



#[cfg(test)]
mod tests {

    
}