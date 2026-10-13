// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, Numeric, TimeFloat};
use common::time_temperature::TimeTemperature;
use io::chrono_inputs::RandomWalkInputs;
use rand::{Rng, RngExt};

const MIN_TIME_STEP: TimeFloat = 1.0;

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

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TimeTempProfile {
    pub times: Vec<TimeFloat>,
    pub temps: Vec<Float>,
}

impl TimeTempProfile {
    pub fn generate_profile(
        &self,
        time_unit: TimeUnit,
        temp_unit: TemperatureUnit,
    ) -> Result<TimeTemperature, String> {
        TimeTemperature::new(self.times.clone(), self.temps.clone(), time_unit, temp_unit)
    }
}

#[derive(Debug, Clone)]
pub struct RandomWalkConstructionParameters {
    pub time_unit: TimeUnit,
    pub temp_unit: TemperatureUnit,
    pub step_max: usize,
    pub min_time: TimeFloat,
    pub min_temp: Float,
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
            min_time: inputs.min_time,
            min_temp: inputs.min_temp,
            time_start_end_range: vec![t_min, t_max],
            temp_start_range: vec![temp_start_min, temp_start_max],
            temp_end_range: vec![temp_end_min, temp_end_max],
            temp_range: vec![temp_min, temp_max],
            monotonic: inputs.monotonic,
        })
    }

    /// Create a profile using the supplied random-number generator.
    pub fn create_new_profile(&self, rng: &mut impl Rng) -> Result<TimeTempProfile, String> {
        let reverse = self.time_unit.is_ka_or_ma();
        let (start_time, end_time) = if reverse {
            (
                sample_positive(
                    self.time_start_end_range[0],
                    self.time_start_end_range[1],
                    rng,
                ),
                0.0,
            )
        } else {
            (
                0.0,
                sample_positive(
                    self.time_start_end_range[0],
                    self.time_start_end_range[1],
                    rng,
                ),
            )
        };

        let start_temp = sample_inclusive(self.temp_start_range[0], self.temp_start_range[1], rng);
        let end_temp = sample_inclusive(self.temp_end_range[0], self.temp_end_range[1], rng);

        let mut times = vec![start_time];
        let mut temps = vec![start_temp];
        let time_direction = (end_time - start_time).signum();
        let temp_direction = (end_temp - start_temp).signum();

        for _ in 0..self.step_max {
            let current_time = *times.last().expect("start time is present");
            let current_temp = *temps.last().expect("start temperature is present");
            let time_remaining = (end_time - current_time).abs();
            let temp_remaining = (end_temp - current_temp).abs();

            // The final endpoint may be closer than the minimum spacing. Stop
            // adding interior points once either endpoint is within that range.
            if time_remaining < 2.0 * self.min_time || temp_remaining < 2.0 * self.min_temp
            {
                break;
            }

            let next_time = self.sample_interior_time(current_time, end_time, time_direction, rng);
            let Some(next_temp) =
                self.sample_interior_temperature(current_temp, end_temp, temp_direction, rng)
            else {
                break;
            };

            times.push(next_time);
            temps.push(next_temp);
        }

        times.push(end_time);
        temps.push(end_temp);

        Ok(TimeTempProfile { times, temps })
    }

    fn sample_interior_temperature(
        &self,
        current: Float,
        end: Float,
        direction: Float,
        rng: &mut impl Rng,
    ) -> Option<Float> {
        if self.monotonic {
            let low = if direction >= 0.0 {
                current + self.min_temp
            } else {
                end +  self.min_temp
            };
            let high = if direction >= 0.0 {
                end -  self.min_temp
            } else {
                current -  self.min_temp
            };
            return (low <= high).then(|| sample_inclusive(low, high, rng));
        }

        // Non-monotonic walks may move in either direction, but every
        // intermediate change must still meet the one-degree minimum.
        for _ in 0..100 {
            let candidate = sample_inclusive(self.temp_range[0], self.temp_range[1], rng);
            if (candidate - current).abs() <  self.min_temp {
                continue;
            }
            if (candidate - end).abs() <=  self.min_temp {
                return None;
            }
            return Some(candidate);
        }
        None
    }


    fn sample_interior_time(
        &self,
        current: TimeFloat,
        end: TimeFloat,
        direction: TimeFloat,
        rng: &mut impl Rng,
    ) -> TimeFloat {
        // Keep at least one configured time unit between the new point and both
        // its predecessor and the final endpoint whenever an interior point fits.
        let (low, high) = if direction >= 0.0 {
            (current + self.min_time, end - self.min_time)
        } else {
            (end + self.min_time, current - self.min_time)
        };
        sample_inclusive(low, high, rng)
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
    pub fn new_random_walk(inputs: &RandomWalkInputs) -> Result<Self, String> {
        Ok(Self::RandomWalk {
            inputs: RandomWalkConstructionParameters::new(inputs)?,
        })
    }

    pub fn create_new_profile(&self, rng: &mut impl Rng) -> Result<TimeTempProfile, String> {
        match self {
            ProfileType::RandomWalk { inputs, .. } => {
                return inputs.create_new_profile(rng);
            }
            ProfileType::RJMCMC { .. } => {
                return Err("Not currently implemented RJMCMC".to_string());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

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
        assert!(profile
            .temps
            .iter()
            .all(|temp| (0.0..=100.0).contains(temp)));
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

        assert!(profile
            .temps
            .iter()
            .all(|temp| (0.0..=100.0).contains(temp)));
    }
}
