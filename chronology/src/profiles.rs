// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Construct, perturb, and validate candidate time-temperature histories.
//!
//! Independent generation samples constraint boxes and recursively subdivides
//! intervals. Walk generation inserts/removes points or perturbs coordinates.
//! Proposal generation is separate from physical simulation and statistical acceptance.



use common::constants::time::TimeUnit;
use common::numeric::{Float, Numeric, TimeFloat};
use common::time_temperature::TimeTemperature;
use common::random::{generate_half_closed_random, generate_half_open_random, generate_normal_random, generate_open_random, generate_weighted_index_random};
use io::chrono_inputs::{ChronologyInputs};

use rand::Rng;

/// Return two numeric bounds in ascending order.
fn ordered_pair<T: Numeric>(first: T, second: T) -> (T, T) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

/// Draw from the closed range, returning its endpoint directly for a zero-width range.
fn sample_inclusive<T: Numeric>(range: (T,T), rng: &mut impl Rng) -> T 
{
    if range.0 == range.1 {
        range.0
    } else {
        T::sample_closed(range.0, range.1, rng)
    }
}

/// Draw strictly between the bounds, except that a zero-width range returns its endpoint.
fn sample_exclusive<T: Numeric>(range: (T,T), rng: &mut impl Rng) -> T {
    if range.0 == range.1 {
        range.0
    } else {
        T::sample_open(range.0, range.1, rng)
    }
}

/// Return whether a value lies within the supplied inclusive, ascending bounds.
///
/// This comparison does not reject NaN; callers requiring finite values
/// must validate finiteness separately.
pub fn check_within_bounds<T: Numeric>(to_check: T, range: (T,T)) -> bool {
    if to_check < range.0 {
        return false; 
    }
    if to_check > range.1 {
        return false;
    }
    true
} 

/// Average `reps` draws from a closed range to favour values near its centre.
///
/// `reps` must be positive and the bounds must be ordered.
fn intermediate_sample_inclusive<T: Numeric>(
    reps: u8, 
    range: (T,T), 
    rng: &mut impl Rng ) -> T 
where 
    T: std::iter::Sum<T> + std::ops::Div<Output = T> + From<u8> 
{
    (0..reps).map(|_| sample_inclusive(range,rng)).sum::<T>() /T::from(reps)
}

/// Average `reps` draws from an open range to favour values near its centre.
///
/// `reps` must be positive; equal bounds use the endpoint without sampling.
fn intermediate_sample_exclusive<T: Numeric>(
    reps: u8, 
    range: (T,T), 
    rng: &mut impl Rng ) -> T 
where 
    T: std::iter::Sum<T> + std::ops::Div<Output = T> + From<u8> 
{
    (0..reps).map(|_| sample_exclusive(range,rng)).sum::<T>() /T::from(reps)
}


/// Recursively fill midpoint nodes between two already populated profile endpoints.
///
/// Time samples lie between neighbours. Modes `I` and `D` restrict temperature
/// samples to neighbour bounds; other modes use the global temperature range.
/// Both slices must contain the supplied indices and `random_rep_number`
/// must be positive. Mutates only nodes strictly between `left` and `right`.
fn fill_intermediate_nodes(
    times: &mut [TimeFloat],
    temps: &mut [Float],
    left: usize,
    right: usize,
    random_rep_number: u8,
    monotonic: &char,
    temp_range: &(Float, Float),
    rng: &mut impl Rng,
) {
    // No intermediate nodes remain
    if right - left <= 1 {
        return;
    }

    let mid = left + (right - left) / 2;

    // Generate midpoint values
    times[mid] = intermediate_sample_exclusive(
        random_rep_number,
        ordered_pair(times[left], times[right]),
        rng,
    );

    temps[mid] = match monotonic{ 
        'D' => {
                intermediate_sample_inclusive(random_rep_number,
                                            ordered_pair(temps[right], temps[left]),
                                            rng)
        },
        'I' => {
                intermediate_sample_inclusive(random_rep_number,
                                      ordered_pair(temps[left], temps[right]),
                                    rng)
        },
        _ => {
              intermediate_sample_inclusive(random_rep_number, *temp_range, rng)
            }
    };
    
    // Fill the left half
    fill_intermediate_nodes(
        times, temps, left, mid, random_rep_number, monotonic, temp_range, rng,);

    // Fill the right half
    fill_intermediate_nodes(
        times, temps, mid, right, random_rep_number, monotonic, temp_range, rng,
    );
}


#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
/// Paired time and temperature control points for a candidate thermal history.
///
/// Vectors must have equal lengths and at least two points before simulation.
/// Times follow the selected unit: increasing elapsed time or decreasing
/// geological age. The record itself stores no units and does not validate on construction.
///
/// # Example
///
/// ```
/// use chronology::profiles::TimeTempProfile;
/// use common::constants::time::TimeUnit;
///
/// let history = TimeTempProfile {
///     times: vec![0.0, 10.0],
///     temps: vec![300.0, 320.0], // Kelvin for the forward model.
/// };
/// assert_eq!(history.k(), 1);
/// assert_eq!(history.age(), 10.0);
/// let model_profile = history.generate_profile(TimeUnit::Second)?;
/// assert_eq!(model_profile.current_temperature(), 300.0);
/// # Ok::<(), String>(())
/// ```
pub struct TimeTempProfile {
    /// Time coordinates in the profile time unit, paired with temperature coordinates.
    pub times: Vec<TimeFloat>,
    /// Temperature coordinates paired with times; the record itself carries no unit metadata.
    pub temps: Vec<Float>,
}
impl io::chrono_plotting::ProfileCoordinates for TimeTempProfile {
    /// Borrow the paired coordinates for plotting without copying the profile.
    fn coordinates(&self) -> io::chrono_plotting::ProfilePath<'_> {
        io::chrono_plotting::ProfilePath {
            times: &self.times,
            temperatures: &self.temps,
        }
    }
}
impl TimeTempProfile {

    /// Clone the control points into a validated, piecewise-linear simulation profile.
    ///
    /// `time_unit` selects forward or geological time ordering. Temperatures
    /// are passed through and must already use the units expected by the model
    /// (Kelvin for physical rate calculations). Invalid profile data returns an error.
    pub fn generate_profile(
        &self,
        time_unit: TimeUnit,
    ) -> Result<TimeTemperature, String> 
    {
        TimeTemperature::new(self.times.clone(), self.temps.clone(), time_unit)
    }
    /// Return the number of segments (`times.len() - 1`); requires a non-empty profile.
    pub fn k(&self) -> usize {
        self.times.len()-1
    }
    /// Return the absolute time span between the first and last control points.
    ///
    /// The result uses the profile time unit. Requires a non-empty time vector.
    pub fn age(&self,) -> TimeFloat{
        let i = self.times.len()-1;
        (self.times[0]-self.times[i]).abs()
    }

}

#[derive(Debug, Clone)]
/// Bounds, constraint boxes, and recursive subdivision rules for independent profiles.
///
/// Constraint endpoints are sampled first; intermediate nodes are recursively
/// inserted between them. Generation and constraint validation are separate steps.
pub struct RandomProfileConstructor{
    /// Allowed temperature bounds in profile units.
    pub temp_range: (Float, Float),
    /// Allowed bounds for every profile time coordinate, in the setup time unit.
    pub age_range: (TimeFloat, TimeFloat),
    /// Exclusive upper bound on the number of profile segments accepted by validation.
    pub step_max: usize,
    /// Strict minimum directed time difference between adjacent nodes.
    pub min_time: TimeFloat,
    /// Strict minimum absolute temperature difference between adjacent nodes.
    pub min_temp: Float,
    /// Total number of nodes generated after subdivision.
    pub profile_length:usize,
    /// Number of interior nodes inserted between consecutive constraint endpoints.
    pub intermediate_node_num:  usize,
    /// Node indices occupied by the sampled box endpoints.
    pub constraint_indices: Vec<usize>,
    /// Constraint rectangles as paired time bounds and temperature bounds.
    pub constraint_boxes: Vec< ((TimeFloat,TimeFloat), (Float,Float))>,
    /// Number of uniform draws averaged for each intermediate coordinate.
    pub random_rep_number: u8, 
    /// Per-segment modes: `I` increasing, `D` decreasing, or `V` unrestricted.
    pub monotonicity: Vec<char>,
}

impl RandomProfileConstructor{

    /// Read constraint boxes and subdivision rules from chronology inputs.
    ///
    /// Adds a fixed time-zero endpoint, placing it last for geological ages.
    /// Returns errors for mismatched monotonicity counts or unsupported mode
    /// characters (`D`, `I`, `V` for monotonicity; `I`, `E`, `G` for sampling).
    pub fn new(inputs: &ChronologyInputs) -> Result<Self, String> {
        if inputs.random_profiles.monotonicity.len() != inputs.random_profiles.constraint_box_num {
            return Err(format!("Number of constraint boxes must be the same las the monotonicity vector length they had {}, {} respecively",inputs.random_profiles.constraint_box_num, inputs.random_profiles.monotonicity.len()) );
        }
        if !inputs.random_profiles.monotonicity.iter().all(|&c| matches!(c, 'D' | 'I' | 'V')) {
            return Err("Invalid character: only D, I, V are allowed".into());
        }
        let random_rep_number:u8 = match inputs.random_profiles.random_segment_kind {
            'I' => 1,
            'E' => 2,
            'G' => 4,
            _ => return Err(format!("Intermediate Random segment has to one of I, E, G but actual value was {} ",inputs.random_profiles.random_segment_kind))

        };

        let intermediate_node_num: usize = (2_u32.pow(inputs.random_profiles.intermediate_node_halving_num)-1) as usize;
        
        let profile_length: usize = (intermediate_node_num + 1) * inputs.random_profiles.constraint_box_num + 1;
        
        let constraint_indices: Vec<usize> = (0..profile_length)
                                            .step_by(intermediate_node_num + 1)
                                            .collect();

        let fixed_point: ((TimeFloat,TimeFloat), (Float,Float)) = ((0.0,0.0), 
                                                                    ordered_pair(inputs.constraints.fixed_temp_point_range.0, 
                                                                                inputs.constraints.fixed_temp_point_range.1) 
                                                                  );
        let mut constraint_boxes: Vec<((TimeFloat,TimeFloat), (Float,Float))> = Vec::with_capacity(profile_length);
        
        if inputs.setup.time_unit.is_ka_or_ma() {
            for ((t_1,t_2),(temp_1, temp_2)) in inputs.random_profiles.constraint_boxes.iter().copied(){
                constraint_boxes.push((ordered_pair(t_1, t_2), (ordered_pair(temp_1, temp_2)))); 
            }
            constraint_boxes.push(fixed_point);
        } else {
            constraint_boxes.push(fixed_point);
            for ((t_1,t_2),(temp_1, temp_2)) in inputs.random_profiles.constraint_boxes.iter().copied(){
                constraint_boxes.push((ordered_pair(t_1, t_2), (ordered_pair(temp_1, temp_2)))); 
            }
        };
       
        Ok( Self{
            temp_range: ordered_pair(inputs.constraints.temp_range.0,inputs.constraints.temp_range.1),
            age_range: ordered_pair(inputs.constraints.age_range.0, inputs.constraints.age_range.1),
            step_max: inputs.constraints.step_max,
            min_time: inputs.constraints.min_time,
            min_temp: inputs.constraints.min_temp,
            profile_length: profile_length,
            intermediate_node_num: intermediate_node_num,
            constraint_indices: constraint_indices,
            constraint_boxes: constraint_boxes,
            random_rep_number: random_rep_number,
            monotonicity: inputs.random_profiles.monotonicity.clone(),

        })
    }

    /// Sample box endpoints and recursively fill intermediate profile coordinates.
    ///
    /// Uses this constructor's subdivision and sampling settings. The returned
    /// profile has not yet passed `check_profile`; invalid constraints may still
    /// produce a proposal that must be rejected by the caller.
    pub fn construct_new_profile(&self, rng: &mut impl Rng)-> Result<TimeTempProfile, String> {

        let mut times: Vec<TimeFloat> = vec![0.0; self.profile_length];// Vec::with_capacity(self.profile_length);
        let mut temps: Vec<Float> = vec![0.0; self.profile_length];// Vec::with_capacity(self.profile_length);

        for (&index, &(time, temp)) in self.constraint_indices
                                                                          .iter()
                                                                          .zip(self.constraint_boxes.iter())
        {
            times[index] = sample_inclusive(time, rng);
            temps[index] = sample_inclusive(temp, rng);

        }

        for (indices, monotonic ) in self.constraint_indices
                                            .windows(2)
                                            .zip(self.monotonicity.iter()) 
        {
            fill_intermediate_nodes(&mut times, &mut temps, 
                                    indices[0], indices[1], 
                                    self.random_rep_number, monotonic, &self.temp_range, rng);                            

        }

        Ok(TimeTempProfile { times, temps })

    }

    /// Check structural, finite-value, global-bound, spacing, and monotonicity constraints.
    ///
    /// `box_check` enables endpoint and intermediate box checks. Intermediate
    /// boxes test the presence of a time and a temperature independently.
    /// `None` uses per-segment monotonicity; `Some` applies one mode globally.
    /// The permitted segment count is strictly less than `step_max`.
    pub fn check_profile(&self, profile: &TimeTempProfile, box_check: bool, monotonicity: Option<char>, ) -> bool {
    
        let count = profile.times.len();

        // Check structure before indexing points.
        if count != profile.temps.len() || count < 2 {
            return false
        }

        let steps = count - 1;
        if steps >= self.step_max {
            return false;
        }

        // Check every point against the global bounds.
        for (&time, &temperature) in
            profile.times.iter().zip(&profile.temps)
        {
            if !time.is_finite()
                || !check_within_bounds(time, self.age_range)
                || !temperature.is_finite()
                || !check_within_bounds(temperature, self.temp_range)
            {
                return false;
            }
        }
        if box_check {
            if !check_within_bounds(profile.times[0], self.constraint_boxes[0].0){
                return false;
            }
            if !check_within_bounds(profile.temps[0], self.constraint_boxes[0].1){
                return false;
            }
            let con_count = self.constraint_boxes.len()-1;
            if !check_within_bounds(profile.times[steps], self.constraint_boxes[con_count].0){
                return false;
            }
            if !check_within_bounds(profile.temps[steps], self.constraint_boxes[con_count].1){
                return false;
            }
            for box_index in 1..con_count {
                let time = profile.times.clone();
                let time_con = time.iter()
                    .any(|t| {check_within_bounds(*t, self.constraint_boxes[box_index].0)});
                let temp = profile.temps.clone();
                let temp_con = temp.iter()
                    .any(|t| {check_within_bounds(*t, self.constraint_boxes[box_index].1)});
                if !time_con || !temp_con{
                    return false;
                }
            }
        }
        let reverse_time = if self.constraint_boxes[0].0.0 == 0.0 {
            false
        }else{
            true
        };
        // Check time direction and minimum changes between adjacent points.
        for index in 0..steps {
            let dt = profile.times[index + 1] - profile.times[index];
            let dtemp = (profile.temps[index + 1] - profile.temps[index]).abs();

            let directed_dt = if reverse_time { -dt } else { dt };
            if !directed_dt.is_finite() 
            || directed_dt <= self.min_time 
            ||!dtemp.is_finite() 
            || dtemp <= self.min_temp {
                return false;
            }
        }
        
        match monotonicity {
            None => {
                for (indices, &mode) in self.constraint_indices
                    .windows(2)
                    .zip(&self.monotonicity)
                {
                    if mode == 'V' {
                        continue;
                    }

                    let start = indices[0];
                    let end = indices[1];
                    let temperatures = &profile.temps[start..=end];
            
                    let valid = match mode {
                        'I' => temperatures.windows(2).all(|pair| pair[1] > pair[0]),
                        'D' => temperatures.windows(2).all(|pair| pair[1] < pair[0]),
                        'V' => true,
                        _ => false,
                    };

                    if !valid {
                        return false;
                    }
                }
            },
            Some(monotonicity) => {
                let temperatures = &profile.temps.clone();
                let valid = match monotonicity {
                        'I' => temperatures.windows(2).all(|pair| pair[1] > pair[0]),
                        'D' => temperatures.windows(2).all(|pair| pair[1] < pair[0]),
                        'V' => true,
                        _ => false,
                    };

                    if !valid {
                        return false;
                    }

            }
        }

        return true;


   
}





}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Kind of profile change in the reversible-jump Monte Carlo walk.
pub enum MoveKind{
    /// Insert a new control point into an existing interval.
    Birth,
    /// Delete an interior control point.
    Death,
    /// Perturb one eligible time coordinate.
    Time,
    /// Perturb one eligible temperature coordinate.
    Temp,
}
impl MoveKind {
    /// Parse `birth`, `death`, `time`, or `temp`; return an error for any other name.
    pub fn from_name(choice: &str) -> Result<Self, String> {
        match choice {
            "birth" => Ok(Self::Birth),
            "death" => Ok(Self::Death),
            "time" => Ok(Self::Time),
            "temp" => Ok(Self::Temp),
            _ => Err(format!("Unknown RJMCMC move type: {choice}")),
        }
    }
    /// Map indices 0, 1, 2, and 3 to birth, death, time, and temperature moves.
    ///
    /// # Panics
    /// Panics for indices outside `0..4`.
    pub fn from_index(index: usize) -> Self {
        match index {
            0 => MoveKind::Birth,
            1 => MoveKind::Death,
            2 => MoveKind::Time,
            3 => MoveKind::Temp,
            _ => unreachable!("expected an in index in 0..4")
        }
    }
    /// Return the stable move-array index: birth 0, death 1, time 2, temperature 3.
    pub fn to_index(self) -> usize {
        match self {
            Self::Birth => 0,
            Self::Death => 1,
            Self::Time => 2,
            Self::Temp => 3,
        }
    }

    /// Return the lowercase configuration/progress name for this move.
    pub fn to_name(self) -> &'static str {
        match self {
            Self::Birth => "birth",
            Self::Death => "death",
            Self::Time => "time",
            Self::Temp => "temp",
        }
    }
}


#[derive(Debug, Clone)]
/// Normalized move weights with alternatives at the profile-size limits.
///
/// Vector order is birth, death, time, temperature. At one segment death is
/// disabled; at `k_max` segments birth is disabled.
pub struct MoveProbabilities {
    /// Move weights in birth, death, time, temperature order.
    pub all: Vec<Float>,
    /// Normalized move weights with birth probability set to zero.
    pub no_birth: Vec<Float>,
    /// Normalized move weights with death probability set to zero.
    pub no_death: Vec<Float>,
    /// Maximum segment count at which further birth moves are disabled.
    pub k_max: usize, 
}
impl MoveProbabilities {

    /// Normalize move weights for ordinary, birth-disabled, and death-disabled proposals.
    ///
    /// Weights must be non-negative and every required sum must be positive.
    /// `k_max` is the maximum segment count used when selecting a move.
    pub fn new(birth: Float, death: Float, time: Float, temp: Float , k_max: usize) -> Self{
        let total = birth + death + time + temp;
        let all:Vec<Float> = vec![ birth/total, death/total, time/total,temp/total];
        let total = death + time + temp;
        let no_birth:Vec<Float> = vec![ 0.0, death/total, time/total,temp/total];
        let total = birth + time + temp;
        let no_death:Vec<Float> = vec![ birth/total, 0.0, time/total,temp/total];
        Self{
            all: all,
            no_birth: no_birth,
            no_death: no_death,
            k_max: k_max
        }
    }
    /// Sample a move, disabling death at one segment and birth at `k_max`.
    ///
    /// Weights must define a valid distribution and the segment count must stay
    /// within the configured limits.
    pub fn choose_step(&self, rng: &mut impl Rng, k:usize) -> MoveKind{

        let index = if k == 1 {
            generate_weighted_index_random(&self.no_death,rng)
        } else if k == self.k_max {
            generate_weighted_index_random(&self.no_birth,rng)
        } else {
            generate_weighted_index_random(&self.all, rng)
        };
        MoveKind::from_index(index)
    }

}

#[derive(Debug, Clone)]
/// Proposal widths in profile time and temperature units.
///
/// Time and temperature perturbations are Gaussian; birth uses scaled uniform
/// draws for interpolation and temperature displacement.
pub struct MoveSigmas{
    /// Time interpolation scale used when inserting a node.
    pub time_birth: TimeFloat,
    /// Temperature displacement width used when inserting a node.
    pub temp_birth: Float,
    /// Gaussian time perturbation width in profile time units.
    pub time: TimeFloat, 
    /// Gaussian temperature perturbation width in profile temperature units.
    pub temp: Float

}
impl MoveSigmas {

    /// Store time/temperature widths for birth and single-node perturbations.
    ///
    /// Widths use the same units as the profile coordinates; no validation is performed.
    pub fn new(time_birth: TimeFloat, temp_birth:Float, time: TimeFloat, temp:Float) -> Self{  
        Self { time_birth, temp_birth, time, temp }
    }

    /// Clone temperatures and add a Gaussian perturbation at index `k`.
    ///
    /// Requires a valid index; bounds are checked later by profile validation.
    pub fn temp(&self, k:usize, current_temps: &Vec<Float>, rng: &mut impl Rng) -> Vec<Float> {
        let mut temps = current_temps.clone();
        temps[k] += self.temp*generate_normal_random(rng);
        return temps;
    }

    /// Clone times and add a Gaussian perturbation at index `k`.
    ///
    /// Requires a valid index; ordering and bounds are checked later.
    pub fn time(&self, k:usize, current_times: &Vec<TimeFloat> , rng: &mut impl Rng) -> Vec<TimeFloat> {
        let mut times = current_times.clone();
        times[k] += self.time*generate_normal_random(rng);
        return times;
    }

    /// Insert a sampled point into the interval adjacent to anchor `k`.
    ///
    /// Returns cloned, extended coordinate vectors. The neighbour is `k + 1`
    /// when the fixed endpoint is last and `k - 1` otherwise. Indices must be
    /// valid; the proposed time and temperature still require constraint checks.
    pub fn birth(&self, k:usize, current_times: &Vec<TimeFloat>, 
        current_temps: &Vec<Float>, fixed_point_last: bool, rng: &mut impl Rng
    ) -> (Vec<TimeFloat>, Vec<Float>) 
    {

            let t_rand: TimeFloat = generate_open_random(0.0,1.0,rng) * self.time_birth;
            let temp_rand:Float = generate_open_random(-0.5, 0.5, rng) * self.temp_birth;

            let (neighbour, insertion_index) = if fixed_point_last {
                (k + 1, k + 1)
            } else {
                (k - 1, k)
            };
            let new_time = current_times[k] + t_rand * (current_times[neighbour] - current_times[k]);
            let new_temp = current_temps[k] + t_rand * (current_temps[neighbour] - current_temps[k]) - temp_rand;
            
            let mut times = current_times.clone();
            let mut temps = current_temps.clone();
            times.insert(insertion_index, new_time);
            temps.insert(insertion_index, new_temp);
            
            (times, temps)
    }
    /// Clone both coordinate vectors and remove node `k` from each.
    ///
    /// The caller must select a valid removable interior node.
    pub fn death(&self, k:usize, current_times: &Vec<TimeFloat>, 
        current_temps: &Vec<Float> )-> (Vec<TimeFloat>, Vec<Float>)  
    {
        let mut times = current_times.clone();
        let mut temps = current_temps.clone();
        times.remove(k);
        temps.remove(k);
        (times, temps)
    }
    
}

#[derive(Debug, Clone)]
/// Constraints and proposal settings for modifying the current thermal history.
pub struct RandomWalkConstructor {
    /// Allowed bounds for every profile time coordinate, in the setup time unit.
    pub age_range: (TimeFloat, TimeFloat),
    /// Allowed temperature bounds in profile units.
    pub temp_range: (Float, Float),
    /// Whether the protected time-zero endpoint is the last node (geological time).
    pub fixed_point_last: bool,
    /// Whether walk proposals must continue to satisfy box constraints.
    pub maintain_constraint_boxes: bool,
    /// Temperature ordering mode: `I` increasing, `D` decreasing, or `V` unrestricted.
    pub monotonicity: char,
    /// Move-selection distributions including profile-size restrictions.
    pub move_probs: MoveProbabilities,
    /// Widths used by birth and coordinate-perturbation proposals.
    pub move_sigmas: MoveSigmas,
}

impl RandomWalkConstructor{

    /// Build walk bounds, move weights, and proposal scales from chronology inputs.
    ///
    /// Geological time units place the fixed time-zero point last. This method
    /// currently wraps construction in `Ok` without validating weights or scales.
    pub fn new(inputs: &ChronologyInputs) -> Result<Self, String> {
        let move_probs = MoveProbabilities::new(
                                                inputs.random_walk.p_birth, 
                                                inputs.random_walk.p_death, 
                                                inputs.random_walk.p_time, 
                                                inputs.random_walk.p_temp, 
                                                inputs.random_walk.k_max);

        let move_sigmas = MoveSigmas::new(
                                        inputs.random_walk.s_time_birth, 
                                        inputs.random_walk.s_temp_birth, 
                                        inputs.random_walk.s_time, 
                                        inputs.random_walk.s_temp);

        Ok(
            Self{
                age_range: ordered_pair(inputs.constraints.age_range.0,
                                         inputs.constraints.age_range.1),
                temp_range: ordered_pair(inputs.constraints.temp_range.0, 
                                        inputs.constraints.temp_range.1),
                fixed_point_last: inputs.setup.time_unit.is_ka_or_ma(),
                maintain_constraint_boxes: inputs.random_walk.maintain_constraint_boxes,
                monotonicity: inputs.random_walk.monotonicity,
                move_probs: move_probs,
                move_sigmas: move_sigmas,
            }
        )
    }
    /// Replace all walk settings, including the proposal scales tuned during burn-in.
    pub fn update_self(&mut self, new: RandomWalkConstructor) {
        self.age_range = new.age_range;
        self.temp_range= new.temp_range;
        self.fixed_point_last = new.fixed_point_last;
        self.maintain_constraint_boxes = new.maintain_constraint_boxes;
        self.monotonicity = new.monotonicity;
        self.move_probs = new.move_probs;
        self.move_sigmas = new.move_sigmas;
    }

    /// Select an eligible node or birth anchor while protecting the fixed endpoint.
    ///
    /// Birth intervals are weighted by their duration, death chooses an interior
    /// node, and perturbations exclude the fixed time-zero endpoint.
    fn pick_k(
        &self, 
        current_profile: &TimeTempProfile, 
        move_type: &MoveKind, 
        rng: &mut impl Rng
    ) -> usize
    {   
        let k = current_profile.k();
        match move_type {
            MoveKind::Birth => {
                let time_gaps: Vec<Float> = current_profile.times
                                                    .windows(2)
                                                    .map(|pair| (pair[1] - pair[0]) .abs())
                                                    .collect();

                let gap = generate_weighted_index_random(&time_gaps,rng);
                // Anchor on the left for right insertion, or on the right for left insertion.
                if self.fixed_point_last { gap } else { gap + 1 }
            },
            MoveKind::Death => { 
                generate_open_random(0, k, rng)
            },
            _ => {
                if self.fixed_point_last{
                    generate_half_open_random(0, k, rng)
                } else {
                    generate_half_closed_random(0, k, rng)
                }
            }
        }



    }

    /// Propose one birth, death, time, or temperature move from the current profile.
    ///
    /// Returns `(trial, move_kind, changed_index)` without accepting or validating
    /// the proposal. For birth, the index identifies the inserted point; for death,
    /// it identifies the removed point in the old profile. The current profile is
    /// left intact and must contain at least two points.
    pub fn construct_new_profile(&self, current_profile: &TimeTempProfile,  rng: &mut impl Rng
    )-> Result<(TimeTempProfile,MoveKind,usize), String> 
    {
        let move_type = self.move_probs.choose_step(rng, current_profile.k());
        let k_choice =  self.pick_k(current_profile, &move_type, rng);
        
        let (times, temps) = match move_type {

            MoveKind::Birth => {
                self.move_sigmas.birth(k_choice, &current_profile.times, &current_profile.temps,
                    self.fixed_point_last, rng)
            }, 
            MoveKind::Death => {
                self.move_sigmas.death(k_choice, &current_profile.times, &current_profile.temps)

            },
            MoveKind::Time => {
                ( self.move_sigmas.time(k_choice, &current_profile.times, rng),
                    current_profile.temps.clone()
                )
            },
            MoveKind::Temp => {
                ( current_profile.times.clone(),
                    self.move_sigmas.temp(k_choice, &current_profile.temps, rng)
                )
            }
        };
        // Acceptance uses the inserted point's index, rather than the birth anchor.
        let changed_index = if move_type == MoveKind::Birth && self.fixed_point_last {
            k_choice + 1
        } else {
            k_choice
        };
        Ok((TimeTempProfile{times,temps},move_type, changed_index))

    }
}


/// Independent profile generator or a walk initialized by an independent generator.
pub enum ProfileType{
    /// Generate each trial independently from constraint boxes.
    Random {
        /// Independent-profile constructor and shared constraint validator.
        constructor: RandomProfileConstructor,
    },
    /// Propose changes to a current profile, using independent generation for initialization.
    Walk   {
        /// Independent-profile constructor and shared constraint validator.
        constructor: RandomProfileConstructor,
        /// Proposal generator that modifies the current accepted profile.
        walker: RandomWalkConstructor,
    }
}

impl ProfileType {
    /// Select independent random generation (`R`) or random-walk generation (`W`).
    ///
    /// Constructs the shared profile constraints and, for walks, proposal settings.
    /// Unsupported generation flags and constructor errors are returned.
    pub fn new(inputs: &ChronologyInputs) -> Result<Self, String> {
        if !matches!(inputs.setup.profile_generation, 'R' | 'W' ) {
            return Err(format!("Invalid profile generation character: only R, and W are allowed, {} was entered",inputs.setup.profile_generation))
        }
        
        let constructor = RandomProfileConstructor::new(inputs)?;
        if inputs.setup.profile_generation == 'R' {
            Ok(ProfileType::Random { constructor })
        } else if inputs.setup.profile_generation == 'W' {
            let walker = RandomWalkConstructor::new(inputs)?;
            Ok(ProfileType::Walk { constructor, walker })
        } else{
            return Err(format!("Invalid profile generation character: only R, and W are allowed, {} was entered",inputs.setup.profile_generation))
        }
    }
    /// Validate a proposal using constraints appropriate to the generation mode.
    ///
    /// Independent profiles use box and per-segment monotonicity checks; walks
    /// use their configured box-preservation flag and global monotonicity.
    pub fn check_profile(&self, profile: &TimeTempProfile) -> bool {
        match self {
            ProfileType::Random { constructor } => {
              
                if constructor.check_profile(&profile, true, None){
                   true
                }else{
                  false
                }
             
            },
            ProfileType::Walk { constructor, walker } => {
                if constructor.check_profile(
                            &profile, 
                            walker.maintain_constraint_boxes,
                            Some(walker.monotonicity))
                {
                    true
                } else {
                    false
                }
            }
        }  
    } 

    /// Generate and validate a proposal, returning `(profile, move_kind, changed_index)`.
    ///
    /// A constraint failure returns `None` for the profile and still consumes a
    /// trial. A walk without a current profile starts from an independent proposal
    /// and returns no move kind. Generation errors propagate.
    pub fn generate_next_profile(&self, current_profile: &Option<TimeTempProfile>, rng: &mut impl Rng
    ) -> Result<(Option<TimeTempProfile>, Option<MoveKind>, usize),String>{

        match self {
            ProfileType::Random { constructor } => {
                let profile = constructor.construct_new_profile(rng)?;
                if constructor.check_profile(&profile, true, None){
                    Ok((Some(profile),None,0))
                }else{
                    Ok((None,None,0))
                }
             
            },
            ProfileType::Walk { constructor, walker } => {
                match current_profile{
                    Some(current_profile) => {
                        let (profile, 
                            kind, 
                            k_choice) = walker.construct_new_profile(current_profile, rng)?;

                        if constructor.check_profile(
                            &profile, 
                            walker.maintain_constraint_boxes,
                            Some(walker.monotonicity)){

                            Ok((Some(profile),Some(kind),k_choice))
                        }else{
                            Ok((None,Some(kind),k_choice))
                        }
                    },
                    None => {
                        let profile = constructor.construct_new_profile(rng)?;
                        if constructor.check_profile(&profile, true, None){
                            Ok((Some(profile),None,0))
                        }else{
                            Ok((None,None,0))
                        } 
                    }
                }
            }
        }

    }

}

////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Verify that walk moves preserve the fixed endpoint in both directions.
    fn walk_moves_preserve_the_fixed_endpoint_in_both_directions() {
        for fixed_point_last in [true, false] {
            let current = TimeTempProfile {
                times: if fixed_point_last { vec![10.0, 4.0, 0.0] } else { vec![0.0, 4.0, 10.0] },
                temps: vec![100.0, 200.0, 300.0],
            };
            for kind in [MoveKind::Birth, MoveKind::Death, MoveKind::Time, MoveKind::Temp] {
                let mut probabilities = vec![0.0; 4];
                probabilities[kind.to_index()] = 1.0;
                let walker = RandomWalkConstructor {
                    age_range: (0.0, 20.0), temp_range: (0.0, 500.0), fixed_point_last,
                    maintain_constraint_boxes: false, monotonicity: 'N',
                    move_probs: MoveProbabilities::new(probabilities[0], probabilities[1],
                        probabilities[2], probabilities[3], 10),
                    move_sigmas: MoveSigmas::new(1.0, 1.0, 1.0, 1.0),
                };
                let mut rng = common::random::get_std_rng_for_rep(0);
                for _ in 0..100 {
                    let (trial, actual_kind, index) = walker.construct_new_profile(&current, &mut rng).unwrap();
                    assert_eq!(actual_kind, kind);
                    let fixed_before = if fixed_point_last { current.times.len() - 1 } else { 0 };
                    let fixed_after = if fixed_point_last { trial.times.len() - 1 } else { 0 };
                    assert_eq!(trial.times[fixed_after], current.times[fixed_before]);
                    assert_eq!(trial.temps[fixed_after], current.temps[fixed_before]);
                    if kind == MoveKind::Birth {
                        assert!(index > 0 && index < trial.times.len() - 1);
                        assert!(trial.times.windows(2).all(|pair|
                            if fixed_point_last { pair[0] > pair[1] } else { pair[0] < pair[1] }));
                        let (times, temps) = walker.move_sigmas.death(index, &trial.times, &trial.temps);
                        assert_eq!(times, current.times);
                        assert_eq!(temps, current.temps);
                    }
                }
            }
        }
    }
}
