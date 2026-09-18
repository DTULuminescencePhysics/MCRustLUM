// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Compact identifiers and constant-time occupancy partitions for physical sites.
//!
//! A [`crate::place_ids::PlaceAvailability`] stores one population as a permutation split into
//! available and unavailable prefixes. Swapping an identifier across that
//! split changes occupancy without scanning the full crystal, which is
//! important when every Monte Carlo event updates multiple site populations.

use crate::numeric::Float;
use serde::{Deserialize, Serialize};
use rand::{Rng, RngExt};

/// Compact zero-based identifier for a trap, hole, or band-tail site.
///
/// The `u16` representation limits each site collection to fewer than 65,536
/// entries and keeps serialized event records small.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct PlaceId(u16);

impl PlaceId {
    /// Convert a zero-based site index into the compact on-disk identifier.
    ///
    /// Returns an error when the index cannot be represented by `u16`.
    pub fn new(index: usize) -> Result<Self, String> {
        let index = u16::try_from(index).map_err(|_| "trap count exceeds u16::MAX".to_string())?;

        Ok(Self(index))
    }

    /// Return this identifier as a zero-based collection index.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Constant-time partition of site identifiers by their current availability.
///
/// Available identifiers occupy `ids[..available_count]`; unavailable ones
/// occupy the remaining suffix. A reverse-position table lets an event move a
/// known identifier across the boundary with a single swap. The physical
/// meaning is population-specific: an available trap is electron occupied,
/// while an available hole site contains a hole that can recombine.
#[derive(Debug)]
pub struct PlaceAvailability {
    /// A permutation containing every [`PlaceId`] exactly once.
    ids: Box<[PlaceId]>,
    /// `positions[position_id]` gives that place's current index in `ids`.
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
    /// Uniformly sample `n` currently unavailable sites without replacement.
    ///
    /// The selected identifiers are moved into the available partition. This
    /// is used to realise fractional initial trap and hole populations.
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
            let selected = rng.random_range(destination..self.ids.len());
            self.swap_positions(destination, selected);
        }

        self.available_count = new_available_end;

        Ok(())
    }

    /// Move every identifier into the unavailable partition.
    ///
    /// For traps this represents an empty occupied-electron set; the meaning
    /// of availability is intentionally supplied by the calling population.
    pub fn mark_all_occupied(&mut self) {
        self.available_count = 0;
    }

    /// Move every identifier into the available partition.
    pub fn mark_all_available(&mut self) {
        self.available_count = self.ids.len();
    }
    /// Return the fixed number of sites represented by this partition.
    pub fn total(&self) -> usize {
        self.ids.len()
    }

    /// Gives Ids available for reaction
    /// i.e. an occupied trap or unoccupied hole
    pub fn available(&self) -> &[PlaceId] {
        &self.ids[..self.available_count]
    }
    /// Gives a vector of indices of the available PlaceIds
    pub fn available_indices_vec(&self) -> Vec<usize> {
        self.available().iter().copied().map(PlaceId::index).collect()
    }

    /// Gives Ids not currently available for reaction
    /// i.e. an unoccupied trap or occupied hole
    pub fn unavailable(&self) -> &[PlaceId] {
        &self.ids[self.available_count..]
    }

    /// Gives a vector of indices of the unavailable PlaceIds
    pub fn unavailable_indices_vec(&self) -> Vec<usize> {
        self.unavailable().iter().copied().map(PlaceId::index).collect()
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
    /// Returns the unavailablity count
    pub fn unavailable_count(&self) -> usize {
        self.total()-self.available_count
    }
    /// Return the available population as a fraction of the total population.
    ///
    /// In the Monte Carlo trap collection, available identifiers are occupied
    /// electron traps, so this value is the simulated trap filling fraction.
    /// An empty collection produces IEEE `NaN` through `0.0 / 0.0`.
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
