// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, Numeric, TimeFloat};
use common::time_temperature::TimeTemperature;
use io::chrono_inputs::RandomWalkInputs;
use rand::distr::Open01;
use rand::{Rng, RngExt};

fn bounds<T: Numeric>(values: &[T], name: &str) -> Result<(T, T), String> {
    let [first, second] = values else {
        return Err(format!("{name} must contain exactly two values"));
    };

    if !first.is_finite() || !second.is_finite() {
        return Err(format!("{name} values must be finite"));
    }

    Ok(ordered_pair(*first, *second))
}

fn positive_bounds<T: Numeric>(values: &[T], name: &str) -> Result<(T, T), String> {
    let (low, high) = bounds(values, name)?;
    if high <= T::zero() {
        return Err(format!("{name} must contain a value greater than zero"));
    }

    // Sampling from an open interval ensures that zero cannot be selected when
    // the user's range crosses zero.
    Ok((low.max(T::zero()), high))
}

fn intersect_bounds<T: Numeric>(
    first: (T, T),
    second: (T, T),
    error: &str,
) -> Result<(T, T), String> {
    let intersection = (first.0.max(second.0), first.1.min(second.1));
    if intersection.0 > intersection.1 {
        Err(error.to_string())
    } else {
        Ok(intersection)
    }
}

fn ordered_pair<T: Numeric>(first: T, second: T) -> (T, T) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn sample_inclusive<T: Numeric>(low: T, high: T, rng: &mut impl Rng) -> T {
    if low == high {
        low
    } else {
        T::random_range(low, high, rng)
    }
}

fn sample_positive<T: Numeric>(low: T, high: T, rng: &mut impl Rng) -> T {
    let zero = T::zero();
    if low > zero {
        return sample_inclusive(low, high, rng);
    }

    // An inclusive range may select zero. Retry, then use the already
    // validated positive upper bound as a deterministic fallback.
    for _ in 0..100 {
        let sample = sample_inclusive(low, high, rng);
        if sample > zero {
            return sample;
        }
    }
    high
}

fn append_random_interior_times(
    times: &mut Vec<TimeFloat>,
    start: TimeFloat,
    end: TimeFloat,
    count: usize,
    rng: &mut impl Rng,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }

    let low = start.min(end);
    let high = start.max(end);
    let width = high - low;
    let mut interior = Vec::with_capacity(count);

    // Open01 excludes both endpoints. Retry the vanishingly unlikely case of
    // duplicate floating-point samples so the ordering remains strictly monotonic.
    for _ in 0..100 {
        interior.clear();
        for _ in 0..count {
            let fraction: TimeFloat = rng.sample(Open01);
            interior.push(low + fraction * width);
        }
        interior.sort_by(TimeFloat::total_cmp);
        if interior.windows(2).all(|pair| pair[0] < pair[1])
            && interior[0] > low
            && interior[count - 1] < high
        {
            if start > end {
                interior.reverse();
            }
            times.extend(interior.iter().copied());
            return Ok(());
        }
    }

    Err("time range is too narrow to generate distinct intermediate times".to_string())
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TimeTempProfile {
    pub times: Vec<TimeFloat>,
    pub temps: Vec<Float>,
}

impl TimeTempProfile {
    pub fn generate_profile(&self, time_unit: TimeUnit, temp_unit: TemperatureUnit) -> Result<TimeTemperature, String> {
        TimeTemperature::new(
            self.times.clone(),
            self.temps.clone(),
            time_unit,
            temp_unit,
        )
    }
}

#[derive(Debug, Clone)]
pub struct RandomWalkConstructionParameters {
    pub time_unit: TimeUnit,
    pub temp_unit: TemperatureUnit,
    pub step_max: usize,
    pub time_start_end_range: Vec<TimeFloat>,
    pub temp_start_range: Vec<Float>,
    pub temp_end_range: Vec<Float>,
    pub temp_range: Vec<Float>,
    pub monotonic: bool,
}

impl RandomWalkConstructionParameters {
    pub fn new(inputs: &RandomWalkInputs) -> Result<Self, String> {
        let reverse = inputs.time_unit.is_ka_or_ma();

        let (t_min, t_max): (TimeFloat, TimeFloat) = if reverse {
            positive_bounds(&inputs.time_start_range, "time_start_range")?
        } else {
            positive_bounds(&inputs.time_end_range, "time_end_range")?
        };
        let (temp_min, temp_max): (Float, Float) = bounds(&inputs.temp_range, "temp_range")?;
        let (temp_start_min, temp_start_max): (Float, Float) = intersect_bounds(
            bounds(&inputs.temp_start_range, "temp_start_range")?,
            (temp_min, temp_max),
            "temp_start_range and temp_range do not overlap",
        )?;
        let (temp_end_min, temp_end_max): (Float, Float) = intersect_bounds(
            bounds(&inputs.temp_end_range, "temp_end_range")?,
            (temp_min, temp_max),
            "temp_end_range and temp_range do not overlap",
        )?;
        Ok(Self {
            time_unit: inputs.time_unit,
            temp_unit: inputs.temp_unit,
            step_max: inputs.step_max,
            time_start_end_range: vec![t_min, t_max],
            temp_start_range: vec![temp_start_min, temp_start_max],
            temp_end_range: vec![temp_end_min, temp_end_max],
            temp_range: vec![temp_min, temp_max],
            monotonic: inputs.monotonic,
        })
    }

    /// Create a profile using the supplied random-number generator.
    pub fn create_new_profile(&self, rng: &mut impl Rng) -> Result<TimeTempProfile, String> {
        let intermediate_steps = rng.random_range(0..=self.step_max);
        let point_count = intermediate_steps + 2;

        let reverse = self.time_unit.is_ka_or_ma();
        let mut times = Vec::with_capacity(point_count);

        if reverse {
            let start = sample_positive(
                self.time_start_end_range[0],
                self.time_start_end_range[1],
                rng,
            );
            times.push(start);
            append_random_interior_times(&mut times, start, 0.0, intermediate_steps, rng)?;
            times.push(0.0);
        } else {
            let end = sample_positive(
                self.time_start_end_range[0],
                self.time_start_end_range[1],
                rng,
            );
            times.push(0.0);
            append_random_interior_times(&mut times, 0.0, end, intermediate_steps, rng)?;
            times.push(end);
        }

        let start_temp = sample_inclusive(self.temp_start_range[0], self.temp_start_range[1], rng);
        let end_temp = sample_inclusive(self.temp_end_range[0], self.temp_end_range[1], rng);
        let mut temps = Vec::with_capacity(point_count);
        temps.push(start_temp);

        if self.monotonic {
            let (low, high) = ordered_pair(start_temp, end_temp);
            for _ in 0..intermediate_steps {
                temps.push(sample_inclusive(low, high, rng));
            }

            if start_temp <= end_temp {
                temps[1..].sort_by(Float::total_cmp);
            } else {
                temps[1..].sort_by(|a, b| b.total_cmp(a));
            }
        } else {
            for _ in 0..intermediate_steps {
                temps.push(sample_inclusive(
                    self.temp_range[0],
                    self.temp_range[1],
                    rng,
                ));
            }
        }
        temps.push(end_temp);

        Ok(TimeTempProfile {
            times,
            temps,
        })
    }
}

pub enum ProfileType {
    RandomWalk {
        inputs: RandomWalkConstructionParameters,
    },
    RJMCMC {
        profile: TimeTempProfile,
    },
}

impl ProfileType {
    pub fn new_random_walk(inputs: &RandomWalkInputs) -> Result<Self,String> {
        Ok(Self::RandomWalk{inputs: RandomWalkConstructionParameters::new(inputs)?})

    }

    pub fn create_new_profile(&self, rng: &mut impl Rng) -> Result<TimeTempProfile, String> {
        match self {
            ProfileType::RandomWalk { inputs, .. } => {return inputs.create_new_profile(rng);},
            ProfileType::RJMCMC { .. } => {return Err("Not currently implemented RJMCMC".to_string());},
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn parameters(time_unit: TimeUnit, monotonic: bool) -> RandomWalkConstructionParameters {
        RandomWalkConstructionParameters {
            time_unit,
            temp_unit: TemperatureUnit::Celsius,
            step_max: 20,
            time_start_end_range: vec![50.0, 100.0],
            temp_start_range: vec![10.0, 20.0],
            temp_end_range: vec![80.0, 90.0],
            temp_range: vec![0.0, 100.0],
            monotonic,
        }
    }

    #[test]
    fn forward_profile_has_strictly_increasing_times() {
        let mut rng = StdRng::seed_from_u64(1);
        let profile = parameters(TimeUnit::Year, false)
            .create_new_profile(&mut rng)
            .unwrap();

        assert_eq!(profile.times[0], 0.0);
        assert!(profile.times.windows(2).all(|pair| pair[0] < pair[1]));
        assert!((50.0..=100.0).contains(profile.times.last().unwrap()));
        assert_eq!(profile.times.len(), profile.temps.len());
    }

    #[test]
    fn reverse_profile_has_strictly_decreasing_times() {
        let mut rng = StdRng::seed_from_u64(2);
        let profile = parameters(TimeUnit::KAnnum, false)
            .create_new_profile(&mut rng)
            .unwrap();

        assert!((50.0..=100.0).contains(&profile.times[0]));
        assert_eq!(*profile.times.last().unwrap(), 0.0);
        assert!(profile.times.windows(2).all(|pair| pair[0] > pair[1]));
    }

    #[test]
    fn monotonic_temperatures_follow_the_endpoint_direction() {
        let mut rng = StdRng::seed_from_u64(3);
        let profile = parameters(TimeUnit::Year, true)
            .create_new_profile(&mut rng)
            .unwrap();

        assert!(profile.temps.windows(2).all(|pair| pair[0] <= pair[1]));
        assert!(
            profile
                .temps
                .iter()
                .all(|temp| (0.0..=100.0).contains(temp))
        );
    }

    #[test]
    fn descending_monotonic_temperatures_are_non_increasing() {
        let mut params = parameters(TimeUnit::Year, true);
        params.temp_start_range = vec![80.0, 90.0];
        params.temp_end_range = vec![10.0, 20.0];
        let mut rng = StdRng::seed_from_u64(4);
        let profile = params.create_new_profile(&mut rng).unwrap();

        assert!(profile.temps.windows(2).all(|pair| pair[0] >= pair[1]));
    }

    #[test]
    fn temperature_endpoint_ranges_are_clipped_to_the_global_range() {
        let mut params = parameters(TimeUnit::Year, false);
        params.temp_start_range = vec![-10.0, 10.0];
        params.temp_end_range = vec![90.0, 110.0];
        params.temp_range = vec![0.0, 100.0];
        let mut rng = StdRng::seed_from_u64(5);
        let profile = params.create_new_profile(&mut rng).unwrap();

        assert!(
            profile
                .temps
                .iter()
                .all(|temp| (0.0..=100.0).contains(temp))
        );
    }
}
