// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Program-wide deterministic random number generation.

use std::sync::{
    LazyLock, Mutex, MutexGuard,
    atomic::{AtomicU64, Ordering},
};

use rand::{SeedableRng, rngs::StdRng};

/// Seed used when [`set_seed`] is not called during program startup.
pub const DEFAULT_SEED: u64 = 0;

static SEED: AtomicU64 = AtomicU64::new(DEFAULT_SEED);
static RNG: LazyLock<Mutex<StdRng>> =
    LazyLock::new(|| Mutex::new(StdRng::seed_from_u64(DEFAULT_SEED)));

/// Set the program-wide random seed and restart the random sequence.
///
/// Call this once near the beginning of the program, before generating any
/// random values. If it is never called, [`DEFAULT_SEED`] is used.
pub fn set_seed(seed: u64) {
    SEED.store(seed, Ordering::Relaxed);
    *rng() = StdRng::seed_from_u64(seed);
}

/// Return the seed most recently supplied to [`set_seed`].
pub fn seed() -> u64 {
    SEED.load(Ordering::Relaxed)
}

/// Lock and return the program-wide random number generator.
///
/// Other modules can use the returned value with methods from [`rand::Rng`].
pub fn rng() -> MutexGuard<'static, StdRng> {
    RNG.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
}
