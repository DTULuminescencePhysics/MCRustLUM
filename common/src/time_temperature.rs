/// This module contains controls the time and temperature progression 

use crate::constants::time;
use crate::numeric::{TimeFloat};
/// Precision used internally by time/temperature profiles.
/// This is intentionally independent of `crate::numeric::Float`, so the rest
/// of the program can use `f32` while time profiles keep enough precision for
/// `ka` and `ma` second values.
/// Time enum can either run forward or backwards in time
#[derive(Debug, Clone, Copy)]
pub enum TimeAdvance {
    Forward {
        start: TimeFloat,
        end: TimeFloat,
        current: TimeFloat,
    },
    Reverse {
        start: TimeFloat,
        end: TimeFloat,
        current: TimeFloat,
    },
}

impl TimeAdvance {
    pub fn forward(end: TimeFloat) -> Self {
        Self::Forward {
            start: 0.0,
            end,
            current: 0.0,
        }
    }

    pub fn reverse(start: TimeFloat) -> Self {
        Self::Reverse {
            start,
            end: 0.0,
            current: start,
        }
    }

    pub fn current(&self) -> TimeFloat {
        match self {
            Self::Forward { current, .. } => *current,
            Self::Reverse { current, .. } => *current,
        }
    }

    pub fn is_forward(&self) -> bool {
        matches!(self, Self::Forward { .. })
    }

    pub fn direction(&self) -> TimeFloat {
        if self.is_forward() {
            1.0
        } else {
            -1.0
        }
    }

    pub fn set_current(&mut self, new_time: TimeFloat) {
        match self {
            Self::Forward { current, .. } => *current = new_time,
            Self::Reverse { current, .. } => *current = new_time,
        }
    }

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

    pub fn is_finished(&self) -> bool {
        match self {
            Self::Forward { current, end, .. } => *current >= *end,
            Self::Reverse { current, end, .. } => *current <= *end,
        }
    }

    pub fn dt_check(&self, dt: &TimeFloat) -> bool {
        match self {
           Self::Forward{ .. }=>  if dt < &0.0 { return false} else{ return true}, 
           Self::Reverse{ .. } => if dt > &0.0 { return false} else{ return true},
        }
       
    }
}

#[derive(Debug, Clone)]
pub struct Step {
    /// dT/dt for each section.
    /// Length is times.len() - 1.
    pub gradients: Vec<TimeFloat>,

    /// Maximum dt for each section such that temperature changes
    /// by no more than 1 degree.
    /// Length is times.len() - 1.
    pub max_dt: Vec<TimeFloat>,

    /// Controls the time advancement
    pub time_advance: TimeAdvance,

    /// Cached maximum allowed dt from the current state.
    /// This is the smaller of:
    /// - max_dt[section_index]
    /// - time until the next section boundary
    pub current_max_dt: TimeFloat,
}
impl Step {

    pub fn new(
        times: &Vec<TimeFloat>,
        temperatures: &Vec<f32>,
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
    
    pub fn dt_check(&self, dt: &TimeFloat) -> bool{
        self.time_advance.dt_check(dt)
    }

    pub fn advance(&mut self, dt: TimeFloat, section_index: &usize) -> f32{

        if !self.dt_check(&dt){ 
            return 0.0;
        }
        
        if !self.time_advance.is_finished() {
            let temp_change = self.gradients[*section_index] * dt ;
            self.time_advance.advance(dt);
            return temp_change as f32;
        } else {
            return 0.0;
        }

    }

}

/// TimeTemperature struct that controls the system temperature
#[derive(Debug, Clone)]
pub struct TimeTemperature {
    /// times and temperatures contain the points where time and temperature change
    /// and have a minimum size of two for the start and end points
    pub times: Vec<TimeFloat>, 
    pub temperatures: Vec<f32>,

    /// Controls the changes in time and temperature
    pub step: Step,

    /// Current temperature
    pub current_temperature: f32,

    /// Current section index.
    /// Section `i` is between:
    /// times[i] and times[i + 1]
    pub section_index: usize,


}

impl TimeTemperature {
    pub fn new(
        times: Vec<f32>,
        temperatures: Vec<f32>,
        unit: &str,
    ) -> Result<Self, String> {
        if times.len() != temperatures.len() {
            return Err("times and temperatures must have the same length".to_string());
        }

        if times.len() < 2 {
            return Err("at least two time/temperature points are required".to_string());
        }
        let (times, time_advance) = Self::convert_time_temperature(times,unit)?;
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

    pub fn convert_time_temperature(
        times: Vec<f32>,
        unit: &str
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
    fn fix_duplicate_times(times: &mut Vec<TimeFloat>, sec: TimeFloat) {
        for i in 1..times.len() {
            if times[i] == times[i - 1] {
                times[i] += &sec;
            }
        }
    }

    pub fn current_time(&self) -> TimeFloat {
        self.step.time_advance.current()
    }

    pub fn current_temperature(&self) -> f32 {
        self.current_temperature
    }

    pub fn current_max_dt(&self) -> TimeFloat {
        self.step.current_max_dt
    }

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

    fn assert_close_f32(actual: f32, expected: f32) {
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
            "s",
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
            "s",
        );

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "at least two time/temperature points are required"
        );
    }

    #[test]
    fn new_rejects_unknown_time_units() {
        let result = TimeTemperature::new(
            vec![0.0, 1.0],
            vec![20.0, 30.0],
            "fortnight",
        );

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "unknown time unit: fortnight");
    }

    #[test]
    fn convert_time_temperature_converts_forward_units_to_seconds() {
        let (times, time_advance) = TimeTemperature::convert_time_temperature(
            vec![0.0, 1.0, 2.0],
            "h",
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
            "s",
        )
        .unwrap();

        assert!(time_advance.is_forward());
        assert_eq!(times, vec![0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn convert_time_temperature_uses_reverse_time_for_ka_profiles() {
        let (times, time_advance) = TimeTemperature::convert_time_temperature(
            vec![2.0, 1.0, 0.0],
            "ka",
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
            "ma",
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
            "s",
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
            "ka",
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
            "s",
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
            "s",
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
            "s",
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
    fn advance_ignores_dt_with_wrong_direction() {
        let mut profile = TimeTemperature::new(
            vec![0.0, 10.0],
            vec![20.0, 30.0],
            "s",
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
            "ka",
        )
        .unwrap();

        let start_time = profile.current_time();
        profile.advance(-1000000.0);

        assert_close(profile.current_time(), start_time  - 1000000.0);
        assert!(profile.current_temperature() > 20.0);
        assert_close_f32(profile.current_temperature(), 20.00031689);
    }
}