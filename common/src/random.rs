// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Program-wide deterministic random number generation.

use std::sync::atomic::{AtomicU64, Ordering};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Exp1, Normal, StandardNormal, weighted::WeightedIndex};
use rand::distr::{Open01, OpenClosed01};
use rand::{Rng, RngExt};

/// Seed used when [`set_seed`] is not called during program startup.
pub const DEFAULT_SEED: u64 = 0;

/// Base seed from which independent repetition streams are derived.
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

/// Derive a deterministic, well-mixed seed for one repetition index.
///
/// SplitMix64 finalization prevents neighbouring repetition numbers from
/// creating visibly correlated seeds while preserving exact reproducibility
/// for a fixed base seed and index.
pub fn seed_for_rep(rep: usize) -> u64 {
    let rep = rep as u64;
    let mut x = seed().wrapping_add(rep);
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

/// Construct the independent standard random-number stream for a repetition.
pub fn get_std_rng_for_rep(rep: usize) -> StdRng  {
    let random_seed = seed_for_rep(rep);
    StdRng::seed_from_u64(random_seed)

}

/// Draw a unit-rate exponentially distributed random value from `rng`.
///
/// The returned sample is suitable for inverse-rate lifetime calculations and
/// is always finite and non-negative for a valid random-number generator.
pub fn generate_exponential_random(rng: &mut impl Rng) -> f64 {
    Exp1.sample(rng)
}

/// Draw a normally distributed random value from `rng` with standard deviation s.
pub fn generate_scaled_normal_random(rng: &mut impl Rng, s: f64) -> f64 {
    let normal = Normal::new(0.0, s)
        .expect("sigma must be finite and positive");

    normal.sample(rng)
}

/// Draw a normally distributed random value from `rng` with standard deviation 1.0.
pub fn generate_normal_random(rng: &mut impl Rng,) -> f64 {
    StandardNormal.sample(rng)
}

/// Draws an index from a weighted distribution in input vector weights.
pub fn generate_weighted_index_random(weights: &[f64], rng: &mut impl Rng,) -> usize {
    let dist = WeightedIndex::new(weights).unwrap();
    dist.sample(rng)
}

/// Open sample (a,b)
pub trait OpenRandom: Sized{
    /// Creates the actual OpenRandom, Sample (a, b).
    fn sample_open(a: Self, b: Self, rng: &mut impl Rng,) -> Self;
}
/// Closed sample [a,b]
pub trait ClosedRandom: Sized{
    /// Creates the actual ClosedRandom Sample [a, b].
    fn sample_closed(a: Self, b: Self, rng: &mut impl Rng,) -> Self;
}
/// Half open on right side [a,b)
pub trait HalfOpenRandom: Sized{
    /// Creates the actual HalfOpenRandom Sample [a, b).
    fn sample_half_open(a: Self, b: Self, rng: &mut impl Rng,) -> Self;
}
/// Hald open on left side (a,b]
pub trait HalfClosedRandom: Sized{
    /// Creates the actual HalfClosedRandom Sample (a, b].
    fn sample_half_closed(a: Self, b: Self, rng: &mut impl Rng,) -> Self;
}
/// Generates open random number
pub fn generate_open_random<T: OpenRandom>(a: T, b: T, rng: &mut impl Rng,) -> T {
    T::sample_open(a, b, rng)
}
/// Generates closed random number
pub fn generate_closed_random<T: ClosedRandom>(a: T, b: T, rng: &mut impl Rng,) -> T {
    T::sample_closed(a, b, rng)
}
/// Generates half open random number
pub fn generate_half_open_random<T: HalfOpenRandom>(a: T, b: T, rng: &mut impl Rng,) -> T {
    T::sample_half_open(a, b, rng)
}
/// Generates half closed random number
pub fn generate_half_closed_random<T: HalfClosedRandom>(a: T, b: T, rng: &mut impl Rng,) -> T {
    T::sample_half_closed(a, b, rng)
}

/// Rules for integers
macro_rules! impl_integer {
    ($($ty:ty),*) => {
        $(
            impl OpenRandom for $ty {
                fn sample_open(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    let start = a.checked_add(1).expect("lower bound overflows");
                    assert!(start < b, "interval contains no integers");
                    rng.random_range(start..b)
                }
            }
            impl ClosedRandom for $ty {
                fn sample_closed(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    assert!(a <= b, "interval is too small or incorrect way round");
                    
                    rng.random_range(a..=b)
                }
            }
            impl HalfOpenRandom for $ty {
                fn sample_half_open(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    rng.random_range(a..b)
                }
            }
            impl HalfClosedRandom for $ty {
                fn sample_half_closed(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    let start = a.checked_add(1).expect("lower bound overflows");
                    rng.random_range(start..=b)
                }
            }

        )*
    };
}
impl_integer!(u8, u16, u32, u64, u128, usize,
              i8, i16, i32, i64, i128);

impl OpenRandom for isize {
    fn sample_open(a: Self, b: Self, rng: &mut impl Rng) -> Self {
        let start = a.checked_add(1).expect("lower bound overflows");
        assert!(start < b, "interval contains no integers");

        rng.random_range((start as i128)..(b as i128)) as isize
    }
}
impl ClosedRandom for isize {
    fn sample_closed(a: Self, b: Self, rng: &mut impl Rng) -> Self {
        assert!(a <= b, "interval is too small or incorrect way round");
        rng.random_range((a as i128)..=(b as i128)) as isize
    }
}
impl HalfOpenRandom for isize {
    fn sample_half_open(a: Self, b: Self, rng: &mut impl Rng) -> Self {
        rng.random_range((a as i128)..(b as i128)) as isize
    }
}
impl HalfClosedRandom for isize {
    fn sample_half_closed(a: Self, b: Self, rng: &mut impl Rng) -> Self {
        let start = a.checked_add(1).expect("lower bound overflows");
        rng.random_range((start as i128)..=(b as i128)) as isize
    }
}
/// Rules for Floats
macro_rules! impl_float {
    ($($ty:ty),*) => {
        $(
            impl OpenRandom for $ty {
                fn sample_open(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    assert!(
                        a.is_finite() && b.is_finite() && a < b,
                        "bounds must be finite and increasing"
                    );

                    loop {
                        let u: Self = rng.sample(Open01);
                        let value = (1.0 - u) * a + u * b;

                        // Floating-point rounding can reach an endpoint.
                        if a < value && value < b {
                            return value;
                        }

                        assert!(
                            a.next_up() < b,
                            "interval contains no representable floats"
                        );
                    }
                }
            }
            impl ClosedRandom for $ty {
                fn sample_closed(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    assert!(a <= b, "interval is too small or incorrect way round");
                    
                    rng.random_range(a..=b)
                }
            }
            impl HalfOpenRandom for $ty {
                fn sample_half_open(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    rng.random_range(a..b)
                }
            }
            impl HalfClosedRandom for $ty {
                fn sample_half_closed(a: Self, b: Self, rng: &mut impl Rng) -> Self {
                    assert!(a.is_finite() && b.is_finite() && a < b);
                    loop {
                        let u: Self = rng.sample(OpenClosed01);
                        let value = (1.0 - u) * a + u * b;

                        if a < value && value <= b {
                            return value;
                        }
                    }
                }
            }
        )*
    };
}

impl_float!(f32, f64);




/// Generates a log of uniform random number in the closed interval (0, 1).
pub fn generate_log_open_uniform_random(rng: &mut impl Rng,) -> f64 {
    let a: f64 = rng.sample(Open01);
     a.ln()
}

