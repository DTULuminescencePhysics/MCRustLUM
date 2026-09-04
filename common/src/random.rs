// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Program-wide deterministic random number generation.

use std::sync::atomic::{AtomicU64, Ordering};
use rand::SeedableRng;
use rand::rngs::StdRng;

/// Seed used when [`set_seed`] is not called during program startup.
pub const DEFAULT_SEED: u64 = 0;

static SEED: AtomicU64 = AtomicU64::new(DEFAULT_SEED);

/// Set the program-wide random seed and restart the random sequence.
///
/// Call this once near the beginning of the program, before generating any
/// random values. If it is never called, [`DEFAULT_SEED`] is used.
pub fn set_seed(seed: u64) {
    SEED.store(seed, Ordering::Relaxed);
}

/// Return the seed most recently supplied to [`set_seed`].
pub fn seed() -> u64 {
    SEED.load(Ordering::Relaxed)
}

pub fn seed_for_rep(rep: usize) -> u64 {
    let rep = rep as u64;
    let mut x = seed().wrapping_add(rep);
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

pub fn get_std_rng_for_rep(rep: usize) -> StdRng  {
    let random_seed = seed_for_rep(rep);
    StdRng::seed_from_u64(random_seed)

}