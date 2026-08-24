//! Piecewise-linear control of simulation time and temperature.
//!
//! A [`TimeTemperature`] profile is defined by matching time and temperature
//! control points. Between adjacent points, temperature changes linearly.
//! [`Step`] caches the gradient for each section and provides a signed maximum
//! timestep that limits a temperature change to one degree and prevents a
//! step from crossing the next control-point boundary.
//!
//! Ordinary units run forward from zero. Geological-age units (`ka` and `ma`)
//! run in reverse from the oldest supplied age toward zero, so their valid
//! timesteps and timestep limits are negative.

use crate::constants::time::TimeUnit;
use crate::constants::temperature::TemperatureUnit;
use crate::constants::{time, temperature};
use crate::numeric::{TimeFloat, Float};

/// Direction and current position of a bounded time interval.
///
/// Time is stored as [`TimeFloat`] independently of the temperature precision.
/// This retains useful resolution after `ka` and `ma` values are converted to
/// their much larger representations in seconds.
#[derive(Debug, Clone, Copy)]
pub enum TimeAdvance {
    /// Advances from `start` to `end` using non-negative timesteps.
    Forward {
        start: TimeFloat,
        end: TimeFloat,
        current: TimeFloat,
    },
    /// Advances from `start` to `end` using non-positive timesteps.
    Reverse {
        start: TimeFloat,
        end: TimeFloat,
        current: TimeFloat,
    },
}

impl TimeAdvance {
    /// Create a forward interval starting at zero and ending at `end`.
    pub fn forward(end: TimeFloat) -> Self {
        Self::Forward {
            start: 0.0,
            end,
            current: 0.0,
        }
    }

    /// Create a reverse interval starting at `start` and ending at zero.
    pub fn reverse(start: TimeFloat) -> Self {
        Self::Reverse {
            start,
            end: 0.0,
            current: start,
        }
    }

    /// Return the current time in seconds.
    pub fn current(&self) -> TimeFloat {
        match self {
            Self::Forward { current, .. } => *current,
            Self::Reverse { current, .. } => *current,
        }
    }

    /// Return whether this interval advances toward increasing time.
    pub fn is_forward(&self) -> bool {
        matches!(self, Self::Forward { .. })
    }

    /// Return `1.0` for forward time or `-1.0` for reverse time.
    ///
    /// This sign is used when constructing direction-aware timestep limits.
    pub fn direction(&self) -> TimeFloat {
        if self.is_forward() {
            1.0
        } else {
            -1.0
        }
    }

    /// Replace the current time without clamping or validating it.
    ///
    /// Callers are responsible for keeping `new_time` inside this interval.
    pub fn set_current(&mut self, new_time: TimeFloat) {
        match self {
            Self::Forward { current, .. } => *current = new_time,
            Self::Reverse { current, .. } => *current = new_time,
        }
    }

    /// Restore the interval to its starting time.
    pub fn reset(&mut self) {
        match self {
            Self::Forward { start, current, .. }
            | Self::Reverse { start, current, .. } => *current = *start,
        }
    }

    /// Move by `dt` and clamp the result to the terminal boundary.
    ///
    /// Forward intervals expect a non-negative value and reverse intervals a
    /// non-positive value. The requirement is checked by a debug assertion;
    /// [`Step::advance`] performs the corresponding runtime check.
    pub fn advance(&mut self, dt: TimeFloat) {
        match self {
            Self::Forward { current, end, .. } => {
                debug_assert!(dt >= 0.0, "Forward TimeAdvance requires positive dt");

                *current = (*current + dt).min(*end);
            }

            Self::Reverse { current, end, .. } => {
                debug_assert!(dt <= 0.0, "Reverse TimeAdvance requires negative dt");

                *current = (*current + dt).max(*end);
            }
        }
    }

    /// Return whether the terminal boundary has been reached.
    pub fn is_finished(&self) -> bool {
        match self {
            Self::Forward { current, end, .. } => *current >= *end,
            Self::Reverse { current, end, .. } => *current <= *end,
        }
    }

    /// Check that a timestep has the correct sign for this interval.
    ///
    /// Zero is valid in either direction.
    pub fn dt_check(&self, dt: &TimeFloat) -> bool {
        match self {
           Self::Forward{ .. }=>  if dt < &0.0 { return false} else{ return true}, 
           Self::Reverse{ .. } => if dt > &0.0 { return false} else{ return true},
        }
       
    }
}

/// Precomputed interpolation data and timestep limits for a profile.
///
/// Section `i` spans control points `i` and `i + 1`. All section vectors
/// therefore contain one fewer element than the profile's control-point
/// vectors.
#[derive(Debug, Clone)]
pub struct Step {
    /// Temperature gradient `dT/dt` for each section.
    /// Length is times.len() - 1.
    pub gradients: Vec<TimeFloat>,

    /// Signed timestep for each section that changes temperature by at most
    /// one degree. Constant-temperature sections use signed infinity.
    /// Length is times.len() - 1.
    pub max_dt: Vec<TimeFloat>,

    /// Direction, bounds, and current position of the profile.
    pub time_advance: TimeAdvance,

    /// Cached signed timestep limit at the current profile position.
    ///
    /// Its magnitude is the smaller of the section's one-degree limit and the
    /// distance to the next section boundary. It is zero when the profile is
    /// finished.
    pub current_max_dt: TimeFloat,
}
impl Step {

    /// Precompute gradients and one-degree timestep limits for every section.
    ///
    /// `times` and `temperatures` must have matching lengths of at least two
    /// and must already be ordered in the direction of `time_advance`.
    pub fn new(
        times: &Vec<TimeFloat>,
        temperatures: &Vec<Float>,
        time_advance: TimeAdvance,
    ) -> Self {
        let direction = time_advance.direction();
        let mut gradients = Vec::with_capacity(times.len() - 1);
        let mut max_dt = Vec::with_capacity(times.len() - 1);

        for i in 0..times.len() - 1 {
            let section_dt = times[i + 1] - times[i];
            let section_dtemp = temperatures[i + 1] - temperatures[i];

            let gradient = section_dtemp as TimeFloat / section_dt;
            gradients.push(gradient);

            let section_max_dt = if gradient.abs() > TimeFloat::EPSILON {
                direction * (1.0 / gradient.abs())
            } else {
                direction * TimeFloat::INFINITY
            };

            max_dt.push(section_max_dt);
        }

        let profile = Self {
            gradients,
            max_dt,
            time_advance,
            current_max_dt: 0.0,
        };
        
        return profile
    }

    /// Recalculate the timestep limit for the current position and section.
    ///
    /// This should be called after advancing or changing sections so a caller
    /// can safely choose its next timestep.
    pub fn update_current_max_dt(&mut self, times: &Vec<TimeFloat>, section_index: &usize) {
        if self.time_advance.is_finished() {
            self.current_max_dt = 0.0;
            return;
        }

        let current_time = self.time_advance.current();
        let next_boundary_time = times[section_index + 1];

        let dt_to_boundary = next_boundary_time - current_time;
        let section_max_dt = self.max_dt[*section_index];

        self.current_max_dt = if section_max_dt.abs() <= dt_to_boundary.abs(){
            section_max_dt
        } else {
            dt_to_boundary
        };   
    }
    
    /// Check whether `dt` has the sign required by the profile direction.
    pub fn dt_check(&self, dt: &TimeFloat) -> bool{
        self.time_advance.dt_check(dt)
    }

    /// Advance time within one section and return its temperature change.
    ///
    /// A timestep with the wrong sign, or a step taken after the profile is
    /// finished, returns zero without changing time. Callers should keep the
    /// magnitude within [`Step::current_max_dt`] so this calculation does not
    /// cross a section boundary.
    pub fn advance(&mut self, dt: TimeFloat, section_index: &usize) -> Float{

        if !self.dt_check(&dt){ 
            return 0.0;
        }
        
        if !self.time_advance.is_finished() {
            let temp_change = self.gradients[*section_index] * dt ;
            self.time_advance.advance(dt);
            return temp_change as Float;
        } else {
            return 0.0;
        }

    }

}

/// Stateful piecewise-linear time and temperature profile.
///
/// The first control point supplies the initial temperature. Advancing updates
/// both time and temperature, snaps exactly reached boundaries to their stored
/// temperature, and then refreshes the next permitted timestep.
#[derive(Debug, Clone)]
pub struct TimeTemperature {
    /// Control-point times converted to seconds in profile order.
    ///
    /// Forward profiles are strictly increasing and reverse profiles strictly
    /// decreasing. There are at least two entries.
    pub times: Vec<TimeFloat>, 
    /// Temperatures corresponding one-to-one with [`Self::times`].
    pub temperatures: Vec<Float>,

    /// Interpolation and time-advancement state.
    pub step: Step,

    /// Current interpolated temperature.
    pub current_temperature: Float,

    /// Current section index.
    /// Section `i` is between:
    /// times[i] and times[i + 1]
    pub section_index: usize,


}

impl TimeTemperature {
    /// Build and validate a time/temperature profile.
    ///
    /// `times` and `temperatures` must have the same length and contain at
    /// least two values. Times are converted from `unit` to seconds. Units
    /// recognised as `ka` or `ma` create a reverse profile; all other known
    /// units create a forward profile.
    ///
    /// Returns an error for an unknown unit, mismatched lengths, insufficient
    /// points, or control points that are not strictly ordered after
    /// conversion.
    pub fn new(
        times: Vec<TimeFloat>,
        temperatures: Vec<Float>,
        time_unit: TimeUnit,
        temp_unit:TemperatureUnit
    ) -> Result<Self, String> {
        if times.len() != temperatures.len() {
            return Err("times and temperatures must have the same length".to_string());
        }

        if times.len() < 2 {
            return Err("at least two time/temperature points are required".to_string());
        }
        let (times, time_advance) = Self::convert_time_temperature(times,time_unit)?;
        let is_forward = time_advance.is_forward();

        if is_forward {
            for i in 1..times.len() {
                if times[i] <= times[i - 1] {
                    return Err(
                        "times must be strictly increasing for Forward TimeAdvance".to_string()
                    );
                }
            }
        } else {
            for i in 1..times.len() {
                if times[i] >= times[i - 1] {
                    return Err(
                        "times must be strictly decreasing for Reverse TimeAdvance".to_string()
                    );
                }
            }
        }
        let temperatures = Self::convert_temperatures(temperatures,temp_unit)?;
        let step = Step::new(&times,&temperatures,time_advance);

        let section_index = 0;
        let current_temperature = temperatures[0];

        let mut profile = Self {
            times,
            temperatures,
            step,
            current_temperature,
            section_index,
        };

        profile.step.update_current_max_dt(&profile.times,&profile.section_index);

        Ok(profile)
    }

    /// Convert control-point times to seconds and select their direction.
    ///
    /// Equal neighbouring times are separated by one second in the direction
    /// of travel. This handles distinct input values that collapse to the same
    /// representable second after conversion.
    pub fn convert_time_temperature(
        times: Vec<TimeFloat>,
        unit: TimeUnit
    ) -> Result<(Vec<TimeFloat>, TimeAdvance), String> {
         let mut times = time::convert_to_seconds(unit,times)
            .ok_or_else(|| format!("unknown time unit: {}", unit))?;

       
        let time_advance = if time::is_ka_or_ma(unit) {
            Self::fix_duplicate_times(&mut times, -1.0);
            let end = *times
                .first()
                .ok_or_else(|| "time vector is empty".to_string())?;
            TimeAdvance::reverse(end)
        } else {
            Self::fix_duplicate_times(&mut times, 1.0);
            let end = *times
                .last()
                .ok_or_else(|| "time vector is empty".to_string())?;

            TimeAdvance::forward(end)
        };

        Ok((times, time_advance))
    }
    pub fn convert_temperatures(temperatures: Vec<Float>, unit: TemperatureUnit) -> Result<Vec<Float>, String> {
        Ok(temperature::convert_to_kelvin(unit, temperatures)
            .ok_or_else(|| "temperature conversion failed".to_string())?)
        
     
    }

    /// Separate adjacent duplicate times by the signed number of seconds.
    fn fix_duplicate_times(times: &mut Vec<TimeFloat>, sec: TimeFloat) {
        for i in 1..times.len() {
            if times[i] == times[i - 1] {
                times[i] += &sec;
            }
        }
    }

    /// Return the current profile time in seconds.
    pub fn current_time(&self) -> TimeFloat {
        self.step.time_advance.current()
    }

    /// Return the current interpolated temperature.
    pub fn current_temperature(&self) -> Float {
        self.current_temperature
    }

    /// Return the signed maximum timestep permitted from the current state.
    ///
    /// The value is positive for forward profiles, negative for reverse
    /// profiles, and zero after the profile finishes.
    pub fn current_max_dt(&self) -> TimeFloat {
        self.step.current_max_dt
    }

    /// Restore the profile to its first control point.
    ///
    /// The control points and precomputed interpolation data are reused. Only
    /// the current time, temperature, section, and timestep limit are reset.
    pub fn reset(&mut self) {
        self.step.time_advance.reset();
        self.current_temperature = self.temperatures[0];
        self.section_index = 0;
        self.step
            .update_current_max_dt(&self.times, &self.section_index);
    }

    /// Advance the profile by one timestep.
    ///
    /// The timestep must have the profile's direction and should not exceed
    /// [`Self::current_max_dt`] in magnitude. Reaching a control point snaps
    /// the temperature to its exact stored value and activates the next
    /// section. Zero and wrong-direction timesteps leave the profile unchanged.
    pub fn advance(&mut self, dt: TimeFloat) {
        if dt == 0.0 {
            return;
        }

        self.current_temperature +=  self.step.advance(dt, &self.section_index);
         
        if (self.current_time() - 
        self.times[self.section_index + 1]).abs() <= TimeFloat::EPSILON{ 
            self.current_temperature = self.temperatures[self.section_index + 1];
                if self.section_index + 1 < self.times.len() - 1 {
                    self.section_index += 1;
                }
        }

        self.step.update_current_max_dt(&self.times,&self.section_index);
    }
   
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: TimeFloat = 1.0e-4;

    fn assert_close(actual: TimeFloat, expected: TimeFloat) {
        assert!(
            (actual - expected).abs() <= EPS,
            "expected {expected}, got {actual}"
        );
    }

    fn assert_close_f32(actual: Float, expected: Float) {
        assert!(
            (actual - expected).abs() <= 1.0e-4,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn time_advance_forward_starts_at_zero_and_clamps_at_end() {
        let mut advance = TimeAdvance::forward(10.0);

        assert!(advance.is_forward());
        assert_close(advance.current(), 0.0);
        assert_close(advance.direction(), 1.0);
        assert!(!advance.is_finished());
        assert!(advance.dt_check(&1.0));
        assert!(!advance.dt_check(&-1.0));

        advance.advance(4.0);
        assert_close(advance.current(), 4.0);

        advance.advance(100.0);
        assert_close(advance.current(), 10.0);
        assert!(advance.is_finished());
    }

    #[test]
    fn time_advance_reverse_starts_at_start_and_clamps_at_zero() {
        let mut advance = TimeAdvance::reverse(10.0);

        assert!(!advance.is_forward());
        assert_close(advance.current(), 10.0);
        assert_close(advance.direction(), -1.0);
        assert!(!advance.is_finished());
        assert!(advance.dt_check(&-1.0));
        assert!(!advance.dt_check(&1.0));

        advance.advance(-4.0);
        assert_close(advance.current(), 6.0);

        advance.advance(-100.0);
        assert_close(advance.current(), 0.0);
        assert!(advance.is_finished());
    }

    #[test]
    fn new_rejects_mismatched_time_and_temperature_lengths() {
        let result = TimeTemperature::new(
            vec![0.0, 1.0],
            vec![20.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        );

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "times and temperatures must have the same length"
        );
    }

    #[test]
    fn new_rejects_profiles_with_fewer_than_two_points() {
        let result = TimeTemperature::new(
            vec![0.0],
            vec![20.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        );

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "at least two time/temperature points are required"
        );
    }

    

    #[test]
    fn convert_time_temperature_converts_forward_units_to_seconds() {
        let (times, time_advance) = TimeTemperature::convert_time_temperature(
            vec![0.0, 1.0, 2.0],
            TimeUnit::Hour,
        )
        .unwrap();

        assert!(time_advance.is_forward());
        assert_close(times[0], 0.0);
        assert_close(times[1], 3600.0);
        assert_close(times[2], 7200.0);
        assert_close(time_advance.current(), 0.0);
    }

    #[test]
    fn convert_time_temperature_adds_one_second_to_duplicate_forward_times() {
        let (times, time_advance) = TimeTemperature::convert_time_temperature(
            vec![0.0, 1.0, 1.0, 3.0],
            TimeUnit::Second
        )
        .unwrap();

        assert!(time_advance.is_forward());
        assert_eq!(times, vec![0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn convert_time_temperature_uses_reverse_time_for_ka_profiles() {
        let (times, time_advance) = TimeTemperature::convert_time_temperature(
            vec![2.0, 1.0, 0.0],
            TimeUnit::KAnnum
        )
        .unwrap();

        assert!(!time_advance.is_forward());
        assert!(times[0] > times[1]);
        assert!(times[1] > times[2]);
        assert_close(time_advance.current(), times[0]);
    }

    #[test]
    fn convert_time_temperature_adds_one_second_to_duplicate_reverse_times() {
        let (times, time_advance) = TimeTemperature::convert_time_temperature(
            vec![2.0, 1.0, 1.0, 0.0],
            TimeUnit::MaAnnum
        )
        .unwrap();

        assert!(!time_advance.is_forward());
        assert!(times[0] > times[1]);
        assert!(times[1] > times[2]);
        assert!(times[2] > times[3]);
    }

    #[test]
    fn new_rejects_non_increasing_forward_profiles_after_conversion() {
        let result = TimeTemperature::new(
            vec![0.0, 2.0, 1.0],
            vec![20.0, 30.0, 40.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        );

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "times must be strictly increasing for Forward TimeAdvance"
        );
    }

    #[test]
    fn new_rejects_non_decreasing_reverse_profiles_after_conversion() {
        let result = TimeTemperature::new(
            vec![0.0, 1.0, 2.0],
            vec![20.0, 30.0, 40.0],
            TimeUnit::KAnnum,
            TemperatureUnit::Kelvin,
        );

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "times must be strictly decreasing for Reverse TimeAdvance"
        );
    }

    #[test]
    fn new_builds_forward_profile_with_expected_initial_state() {
        let profile = TimeTemperature::new(
            vec![0.0, 10.0, 20.0],
            vec![20.0, 30.0, 50.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        )
        .unwrap();

        assert_eq!(profile.times, vec![0.0, 10.0, 20.0]);
        assert_eq!(profile.temperatures, vec![20.0, 30.0, 50.0]);
        assert_close(profile.current_time(), 0.0);
        assert_close_f32(profile.current_temperature(), 20.0);
        assert_eq!(profile.section_index, 0);
        assert!(profile.step.time_advance.is_forward());
    }

    #[test]
    fn step_new_calculates_gradients_and_max_dt_for_forward_profile() {
        let times = vec![0.0, 10.0, 20.0];
        let temperatures = vec![20.0, 30.0, 50.0];
        let step = Step::new(&times, &temperatures, TimeAdvance::forward(20.0));

        assert_eq!(step.gradients.len(), 2);
        assert_eq!(step.max_dt.len(), 2);
        assert_close(step.gradients[0], 1.0);
        assert_close(step.gradients[1], 2.0);
        assert_close(step.max_dt[0], 1.0);
        assert_close(step.max_dt[1], 0.5);
    }

    #[test]
    fn update_current_max_dt_uses_smaller_of_temperature_limit_and_boundary_limit() {
        let mut profile = TimeTemperature::new(
            vec![0.0, 0.5, 10.0],
            vec![20.0, 21.0, 30.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        )
        .unwrap();

        profile.step.update_current_max_dt(&profile.times, &profile.section_index);
        assert_close(profile.current_max_dt(), 0.5);
    }

    #[test]
    fn advance_updates_time_temperature_section_and_current_max_dt() {
        let mut profile = TimeTemperature::new(
            vec![0.0, 10.0, 20.0],
            vec![20.0, 30.0, 50.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        )
        .unwrap();

        profile.advance(5.0);
        assert_close(profile.current_time(), 5.0);
        assert_close_f32(profile.current_temperature(), 25.0);
        assert_eq!(profile.section_index, 0);

        profile.advance(5.0);
        assert_close(profile.current_time(), 10.0);
        assert_close_f32(profile.current_temperature(), 30.0);
        assert_eq!(profile.section_index, 1);

        profile.advance(2.0);
        assert_close(profile.current_time(), 12.0);
        assert_close_f32(profile.current_temperature(), 34.0);
        assert_eq!(profile.section_index, 1);
    }

    #[test]
    fn reset_restores_forward_profile_starting_state() {
        let mut profile = TimeTemperature::new(
            vec![0.0, 10.0, 20.0],
            vec![20.0, 30.0, 50.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        )
        .unwrap();

        profile.advance(10.0);
        profile.advance(2.0);
        profile.reset();

        assert_close(profile.current_time(), 0.0);
        assert_close_f32(profile.current_temperature(), 20.0);
        assert_eq!(profile.section_index, 0);
        assert_close(profile.current_max_dt(), 1.0);
    }

    #[test]
    fn advance_ignores_dt_with_wrong_direction() {
        let mut profile = TimeTemperature::new(
            vec![0.0, 10.0],
            vec![20.0, 30.0],
            TimeUnit::Second,
            TemperatureUnit::Kelvin,
        )
        .unwrap();

        profile.advance(-5.0);
        assert_close(profile.current_time(), 0.0);
        assert_close_f32(profile.current_temperature(), 20.0);
    }

    #[test]
    fn reverse_profile_advances_backwards_in_time() {
        let mut profile = TimeTemperature::new(
            vec![2.0, 1.0, 0.0],
            vec![20.0, 30.0, 40.0],
            TimeUnit::KAnnum,
            TemperatureUnit::Kelvin,
        )
        .unwrap();

        let start_time = profile.current_time();
        profile.advance(-1000000.0);

        assert_close(profile.current_time(), start_time  - 1000000.0);
        assert!(profile.current_temperature() > 20.0);
        assert_close_f32(profile.current_temperature(), 20.00031689);
    }

    #[test]
    fn reset_restores_reverse_profile_starting_state() {
        let mut profile = TimeTemperature::new(
            vec![2.0, 1.0, 0.0],
            vec![20.0, 30.0, 40.0],
            TimeUnit::KAnnum,
            TemperatureUnit::Kelvin,
        )
        .unwrap();
        let starting_time = profile.current_time();
        let starting_max_dt = profile.current_max_dt();

        profile.advance(-1000000.0);
        profile.reset();

        assert_close(profile.current_time(), starting_time);
        assert_close_f32(profile.current_temperature(), 20.0);
        assert_eq!(profile.section_index, 0);
        assert_close(profile.current_max_dt(), starting_max_dt);
    }
}
