// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only


use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, Numeric, TimeFloat};
use common::time_temperature::TimeTemperature;
use io::chrono_inputs::{RandomWalkInputs, RJMCMCInputs};
use rand::{Rng, RngExt};
use common::random::{generate_normal_random, generate_weighted_index_random, 
    generate_uniform_one_random, generate_uniform_between_values};


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
pub fn check_within_bounds<T: Numeric>(to_check: T, low: T, high: T) -> bool {
    if to_check < low {
        return false; 
    }
    if to_check > high {
        return false;
    }
    true
} 
pub fn check_if_monotonic<T: Numeric>(to_check: T, before: T, after: T, decreasing: bool) -> bool {

    if decreasing{
        return check_within_bounds(to_check, after, before);
     
    }else {
        return check_within_bounds(to_check, before, after);
    }
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
    ) -> Result<TimeTemperature, String> 
    {
        TimeTemperature::new(self.times.clone(), self.temps.clone(), time_unit, temp_unit)
    }

    pub fn new_with_temperature(
        &self,
        i: usize,
        new_temp: Float,
    ) -> Option<Self> {
        if self.times.len() != self.temps.len() {
            return None
        }
        if i >= self.temps.len() {
            return None
        }

        let mut new_profile = self.clone();
        new_profile.temps[i] = new_temp;
        Some(new_profile)
    }

    pub fn new_with_time(
        &self,
        i: usize,
        new_time: TimeFloat,
    ) -> Option<Self> {
        if self.times.len() != self.temps.len() {
            return None;
        }
        if i >= self.times.len() {
            return None;
        }

        let mut new_profile = self.clone();
        new_profile.times[i] = new_time;
        Some(new_profile)
    }

    pub fn new_with_death(
        &self, i: usize
    ) -> Option<Self> 
    {
        if self.times.len() != self.temps.len() {
            return None;
        }
        if i >= self.times.len() {
            return None;
        }

        let mut new_profile = self.clone();
        new_profile.times.remove(i);
        new_profile.temps.remove(i);
        Some(new_profile)
    }

    pub fn new_with_birth(
        &self,
        i: usize,
        new_time: TimeFloat,
        new_temp: Float,
    ) -> Option<Self> {
        if self.times.len() != self.temps.len() {
            return None;
        }
        if i > self.times.len() {
            return None;
        }

        let mut new_profile = self.clone();
        new_profile.times.insert(i+1, new_time);
        new_profile.temps.insert(i+1, new_temp);
        Some(new_profile)
    }


    pub fn log_prior(
        &self, 
        k_max: usize, 
        temp_min: Float, temp_max:Float, 
        temp_start_min: Float, temp_start_max: Float,
        temp_end_min: Float, temp_end_max: Float 
    ) -> Float 
    {
        let k = self.times.len();
        let (t_min, t_max) = if self.times[0] < self.times[k-1]{
            (self.times[0], self.times[k-1])
        } else{
            (self.times[k-1], self.times[0])
        };

        let mut log_p: Float = - (
                            (k_max as Float).ln() + 
                            ((k as Float)*(t_max-t_min).ln() as Float) +
                            (((k-2) as Float)*(temp_max - temp_min).ln()) +
                            (temp_start_max-temp_start_min).ln() +
                            (temp_end_max-temp_end_min).ln()
                        );

        for i in 0..k {
            log_p +=  ((self.times[i]-self.times[i+1]).abs()).ln() as Float +
                      ((i + 1) as Float).ln();
        }    

        return log_p;

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


#[derive(Debug, Clone)]
pub struct RJMCMCSigmas{
    pub time_birth: TimeFloat,
    pub temp_birth: Float,
    pub time: TimeFloat, 
    pub temp: Float

}
impl RJMCMCSigmas {

    pub fn new(time_birth: TimeFloat, temp_birth:Float, time: TimeFloat, temp:Float) -> Self{  
        Self { time_birth, temp_birth, time, temp }
    }

    pub fn temp(&self, current_temp: Float, rng: &mut impl Rng) -> Float {
        current_temp + self.temp*generate_normal_random(rng)
    }

    pub fn time(&self, current_time: TimeFloat, rng: &mut impl Rng) -> Float {
        current_time + self.time*generate_normal_random(rng)
    }

    pub fn birth(
        &self, 
        before_time: TimeFloat, 
        after_time: TimeFloat, 
        before_temp: Float, 
        after_temp: Float, 
        rng: &mut impl Rng) -> (TimeFloat, Float, Float) {

            let t_rand = generate_uniform_one_random(rng);
            let temp_rand = generate_uniform_between_values(-0.5,0.5, rng);

            (
                before_time - self.time_birth*t_rand*(before_time-after_time),
                before_temp - (before_temp-after_temp)*self.time_birth*t_rand - self.temp_birth*temp_rand,
                t_rand*temp_rand
            )

    }

    pub fn death(&self, rng: &mut impl Rng) -> Float {
        generate_uniform_one_random(rng)*generate_uniform_between_values(-0.5,0.5, rng)
    }
    
}


#[derive(Debug, Clone)]
pub struct RJMCMCConstructionParameters {
    pub time_unit: TimeUnit,
    pub temp_unit: TemperatureUnit,
    pub step_max: usize,
    pub min_time: TimeFloat,
    pub min_temp: Float,
    pub time_start_end_range: Vec<TimeFloat>,
    pub temp_start_range: Vec<Float>,
    pub temp_end_range: Vec<Float>,
    pub temp_range: Vec<Float>,
    pub sigmas: RJMCMCSigmas,
    pub monotonic: bool,
    pub decreasing: bool
}

impl RJMCMCConstructionParameters{

    pub fn new(inputs: &RJMCMCInputs) -> Result<Self, String> {
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
        let sigmas = RJMCMCSigmas::new(inputs.s_time_birth, inputs.s_temp_birth, inputs.s_time, inputs.s_temp);
        let decreasing = if temp_start_min > temp_end_max {
            true
        } else {
            false
        };
        Ok( Self {
            time_unit: inputs.time_unit,
            temp_unit: inputs.temp_unit,
            step_max: inputs.step_max,
            min_time: inputs.min_time,
            min_temp: inputs.min_temp,
            time_start_end_range: vec![t_min, t_max],
            temp_start_range: vec![temp_start_min, temp_start_max],
            temp_end_range: vec![temp_end_min, temp_end_max],
            temp_range: vec![temp_min, temp_max],
            sigmas,
            monotonic: inputs.monotonic,
            decreasing,
        })

    }

    pub fn propose_trial_profile(
        &self, 
        current_profile: &TimeTempProfile, 
        move_choice: &str, 
        rng: &mut impl Rng
    ) -> Result<(Option<TimeTempProfile>, Float), String> 
    {
        
        match move_choice {
            "birth" => Ok(self.alter_birth(current_profile, rng)),
            "death" => Ok(self.alter_death(current_profile, rng)),
            "time"  => Ok(self.alter_time(current_profile, rng)),
            "temp"  => Ok(self.alter_temperature(current_profile, rng)),
            _       => return Err(format!("This outcome should never occur the chosen move is: {move_choice}"))
        }

    }

    pub fn alter_temperature(
        &self,
        current_profile: &TimeTempProfile, 
        rng: &mut impl Rng,
    ) -> (Option<TimeTempProfile>, Float) 
    {
        let k = current_profile.times.len();
        let i: usize = rng.random_range(0..k);

        let new_temp = self.sigmas.temp(current_profile.temps[i], rng);

        let mut lower = self.temp_range[0];
        let mut upper = self.temp_range[1];

        if i == 0 {
            lower = lower.max(self.temp_start_range[0]);
            upper = upper.min(self.temp_start_range[1]);
        }

        if i == k - 1 {
            lower = lower.max(self.temp_end_range[0]);
            upper = upper.min(self.temp_end_range[1]);
        }

        if self.monotonic && k > 1 {
            if i == 0 {
                let next = current_profile.temps[i + 1];
                if self.decreasing {
                    lower = lower.max(next);
                } else {
                    upper = upper.min(next);
                }
            } else if i == k - 1 {
                let prev = current_profile.temps[i - 1];
                if self.decreasing {
                    upper = upper.min(prev);
                } else {
                    lower = lower.max(prev);
                }
            } else {
                let prev = current_profile.temps[i - 1];
                let next = current_profile.temps[i + 1];
                if self.decreasing {
                    lower = lower.max(next);
                    upper = upper.min(prev);
                } else {
                    lower = lower.max(prev);
                    upper = upper.min(next);
                }
            }
        }
        
        if lower <= upper && check_within_bounds(new_temp, lower, upper) {
            (current_profile.new_with_temperature(i, new_temp), 0.0)
        } else {
            (None, 0.0) 
        }

    }

    pub fn alter_time(
        &self,
        current_profile: &TimeTempProfile,
        rng: &mut impl Rng,
    ) -> (Option<TimeTempProfile>, Float) 
    {
        let k = current_profile.times.len();
        let i: usize = if self.time_unit.is_ka_or_ma(){
            rng.random_range(0..k-1)
        } else {
            rng.random_range(1..k)
        };
        
        
        let new_time = self.sigmas.time(current_profile.times[i], rng);
        

        let (lower, upper) = if i == 0 {
            (current_profile.times[i+1],
            self.time_start_end_range[0].max(self.time_start_end_range[1]))
       
        } else if i == k - 1 {
            (self.time_start_end_range[0].min(self.time_start_end_range[1]),
            current_profile.times[i-1])  
        } else {
            if self.time_unit.is_ka_or_ma(){
                (current_profile.times[i+1],current_profile.times[i-1])
            }else{
                (current_profile.times[i-1],current_profile.times[i+1])
            }
        };

        if lower <= upper && check_within_bounds(new_time, lower, upper) {
            let prop = current_profile.new_with_time(i, new_time);
            match prop {
                None => (None, 0.0), 
                Some(profile) => {
                    let k_float = k as TimeFloat - 1.0;
                    let log_prior_ratio = if i == 0 {
                        k_float*(current_profile.times[0].max(current_profile.times[k-1])).ln() -
                        k_float*(profile.times[0].max(profile.times[k-1])).ln() +
                        (profile.times[i]-profile.times[i+1]).abs().ln() -
                        (current_profile.times[i]- current_profile.times[i+1]).abs().ln()
                    } else if i == k - 1 {
                        k_float*(current_profile.times[0].max(current_profile.times[k-1])).ln() -
                        k_float*(profile.times[0].max(profile.times[k-1])).ln() +
                        (profile.times[i-1]-profile.times[i]).abs().ln() + 
                        (current_profile.times[i-1]- current_profile.times[i]).abs().ln() 
                    } else {
                
                    k_float*(current_profile.times[0].max(current_profile.times[k-1])).ln() -
                    k_float*(profile.times[0].max(profile.times[k-1])).ln() +
                    (profile.times[i-1]-profile.times[i]).abs().ln() + 
                    (profile.times[i]-profile.times[i+1]).abs().ln() -
                    (current_profile.times[i-1]- current_profile.times[i]).abs().ln() - 
                    (current_profile.times[i]- current_profile.times[i+1]).abs().ln()
                    };
                    (Some(profile), log_prior_ratio)
                } 
            }
        } else {
            (None, 0.0) 
        }

    }

    pub fn alter_birth(
        &self,
        current_profile: &TimeTempProfile,
        rng: &mut impl Rng,
    ) -> (Option<TimeTempProfile>, Float) 
    {
        let k = current_profile.times.len();
        let time_gaps: Vec<Float> = current_profile.times
                                                    .windows(2)
                                                    .map(|pair| (pair[1] - pair[0]) .abs())
                                                    .collect();

        let i = generate_weighted_index_random(&time_gaps,rng);
        
        let (new_time, new_temp, jacobian) = self.sigmas.birth(
                                                            current_profile.times[i], 
                                                            current_profile.times[i+1], 
                                                            current_profile.temps[i], 
                                                            current_profile.temps[i+1], 
                                                            rng
                                                        );
        
        let (lower, upper) = if self.time_unit.is_ka_or_ma(){
                                            (current_profile.times[i+1],current_profile.times[i])
                                        }else{
                                            (current_profile.times[i],current_profile.times[i+1])
                                        };
        let mut lower_temp = self.temp_range[0];
        let mut upper_temp = self.temp_range[1];
        if self.monotonic && k > 1 {
            let prev = current_profile.temps[i];
            let next = current_profile.temps[i + 1];
                if self.decreasing {
                    lower_temp = lower_temp.max(next);
                    upper_temp = upper_temp.min(prev);
                } else {
                    lower_temp = lower_temp.max(prev);
                    upper_temp = upper_temp.min(next);
                }
           
        }
        if lower <= upper && check_within_bounds(new_time, lower, upper) && 
         lower_temp <= upper_temp && check_within_bounds(new_temp, lower_temp, upper_temp){
            (
                current_profile.new_with_birth(i, new_time, new_temp), 
                (new_time-lower).abs().ln() + (upper - new_time).abs().ln() - 
                (upper-lower).ln() - (self.temp_range[1]-self.temp_range[0]).abs().ln()
                + jacobian.ln() )        
            
        } else {
            (None, 0.0) 
        }

    }

    pub fn alter_death(
        &self,
        current_profile: &TimeTempProfile,
        rng: &mut impl Rng,
    ) -> (Option<TimeTempProfile>, Float) 
    {   
        let k = current_profile.times.len();
        let i: usize = rng.random_range(1..k-1);
       
        let jacobian = self.sigmas.death(rng);
        let old_time = current_profile.times[i];
        ( current_profile.new_with_death(i), 

        (current_profile.times[i-1]-current_profile.times[i+1]).abs().ln() +
        (self.temp_range[1]-self.temp_range[0]).ln() -
        (old_time-current_profile.times[i+1]).abs().ln() -
        (current_profile.times[i-1]-old_time).abs().ln() -jacobian.ln())
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
        inputs: RJMCMCConstructionParameters,
    },
}

impl ProfileType {
    pub fn new_random_walk(inputs: &RandomWalkInputs) -> Result<Self, String> {
        Ok(Self::RandomWalk {
            inputs: RandomWalkConstructionParameters::new(inputs)?,
        })
    }
    pub fn new_rjmcmc(inputs: &RJMCMCInputs) -> Result<Self, String> {
         Ok(Self::RJMCMC {
            inputs: RJMCMCConstructionParameters::new(inputs)?,
        })
    }

    pub fn create_new_profile(&self, rng: &mut impl Rng) -> Result<TimeTempProfile, String> {
        match self {
            ProfileType::RandomWalk { inputs } => {
                return inputs.create_new_profile(rng);
            }
            ProfileType::RJMCMC { inputs } => {
                return inputs.create_new_profile(rng);
            }
        }
    }

    pub fn create_new_trial_profile(
        &self,
        current_profile: &TimeTempProfile, 
        move_choice: &str,  
        rng: &mut impl Rng 
    ) -> Result<(Option<TimeTempProfile>, Float) , String >
    {
        match self {
            ProfileType::RandomWalk { .. } => {
                return Err("Not currently implemented RJMCMC".to_string());
            }
            ProfileType::RJMCMC { inputs, .. } => {
                return inputs.propose_trial_profile(current_profile, move_choice, rng)
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
            min_time: 1.0,
            min_temp: 0.01,
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
