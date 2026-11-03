// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only



use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, Numeric, TimeFloat};
use common::time_temperature::TimeTemperature;
use common::random::{generate_weighted_index_random, 
    generate_normal_random, generate_closed_random,
    generate_open_random};
use io::chrono_inputs::{ChronologyInputs};

use rand::Rng;

fn ordered_pair<T: Numeric>(first: T, second: T) -> (T, T) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn sample_inclusive<T: Numeric>(range: (T,T), rng: &mut impl Rng) -> T 
{
    if range.0 == range.1 {
        range.0
    } else {
        T::sample_closed(range.0, range.1, rng)
    }
}

fn sample_exclusive<T: Numeric>(range: (T,T), rng: &mut impl Rng) -> T {
    if range.0 == range.1 {
        range.0
    } else {
        T::sample_open(range.0, range.1, rng)
    }
}

pub fn check_within_bounds<T: Numeric>(to_check: T, range: (T,T)) -> bool {
    if to_check < range.0 {
        return false; 
    }
    if to_check > range.1 {
        return false;
    }
    true
} 

fn intermediate_sample_inclusive<T: Numeric>(
    reps: u8, 
    range: (T,T), 
    rng: &mut impl Rng ) -> T 
where 
    T: std::iter::Sum<T> + std::ops::Div<Output = T> + From<u8> 
{
    (0..reps).map(|_| sample_inclusive(range,rng)).sum::<T>() /T::from(reps)
}

fn intermediate_sample_exclusive<T: Numeric>(
    reps: u8, 
    range: (T,T), 
    rng: &mut impl Rng ) -> T 
where 
    T: std::iter::Sum<T> + std::ops::Div<Output = T> + From<u8> 
{
    (0..reps).map(|_| sample_exclusive(range,rng)).sum::<T>() /T::from(reps)
}


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
pub struct TimeTempProfile {
    pub times: Vec<TimeFloat>,
    pub temps: Vec<Float>,
}
impl io::chrono_plotting::ProfileCoordinates for TimeTempProfile {
    fn coordinates(&self) -> io::chrono_plotting::ProfilePath<'_> {
        io::chrono_plotting::ProfilePath {
            times: &self.times,
            temperatures: &self.temps,
        }
    }
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
    pub fn k(&self) -> usize {
        self.times.len()-1
    }
    pub fn age(&self,) -> TimeFloat{
        let i = self.times.len()-1;
        (self.times[0]-self.times[i]).abs()
    }

}

#[derive(Debug, Clone)]
pub struct RandomProfileConstructor{
    pub temp_range: (Float, Float),
    pub age_range: (TimeFloat, TimeFloat),
    pub step_max: usize,
    pub min_time: TimeFloat,
    pub min_temp: Float,
    pub profile_length:usize,
    pub intermediate_node_num:  usize,
    pub constraint_indices: Vec<usize>,
    pub constraint_boxes: Vec< ((TimeFloat,TimeFloat), (Float,Float))>,
    pub random_rep_number: u8, 
    pub monotonicity: Vec<char>,
}

impl RandomProfileConstructor{

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

        // Check time direction and minimum changes between adjacent points.
        for index in 0..steps {
            let dt = (profile.times[index + 1] - profile.times[index]).abs();
            let dtemp = (profile.temps[index + 1] - profile.temps[index]).abs();

            if !dt.is_finite() 
            || dt <= self.min_time 
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
pub enum MoveKind{
    Birth,
    Death,
    Time,
    Temp,
}
impl MoveKind {
    pub fn from_name(choice: &str) -> Result<Self, String> {
        match choice {
            "birth" => Ok(Self::Birth),
            "death" => Ok(Self::Death),
            "time" => Ok(Self::Time),
            "temp" => Ok(Self::Temp),
            _ => Err(format!("Unknown RJMCMC move type: {choice}")),
        }
    }
    pub fn from_index(index: usize) -> Self {
        match index {
            0 => MoveKind::Birth,
            1 => MoveKind::Death,
            2 => MoveKind::Time,
            3 => MoveKind::Temp,
            _ => unreachable!("expected an in index in 0..4")
        }
    }
    pub fn to_index(self) -> usize {
        match self {
            Self::Birth => 0,
            Self::Death => 1,
            Self::Time => 2,
            Self::Temp => 3,
        }
    }

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
pub struct MoveProbabilities {
    pub all: Vec<Float>,
    pub no_birth: Vec<Float>,
    pub no_death: Vec<Float>,
    pub k_max: usize, 
}
impl MoveProbabilities {

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
pub struct MoveSigmas{
    pub time_birth: TimeFloat,
    pub temp_birth: Float,
    pub time: TimeFloat, 
    pub temp: Float

}
impl MoveSigmas {

    pub fn new(time_birth: TimeFloat, temp_birth:Float, time: TimeFloat, temp:Float) -> Self{  
        Self { time_birth, temp_birth, time, temp }
    }

    pub fn temp(&self, k:usize, current_temps: &Vec<Float>, rng: &mut impl Rng) -> Vec<Float> {
        let mut temps = current_temps.clone();
        temps[k] += self.temp*generate_normal_random(rng);
        return temps;
    }

    pub fn time(&self, k:usize, current_times: &Vec<TimeFloat> , rng: &mut impl Rng) -> Vec<TimeFloat> {
        let mut times = current_times.clone();
        times[k] += self.time*generate_normal_random(rng);
        return times;
    }

    pub fn birth(&self, k:usize, current_times: &Vec<TimeFloat>, 
        current_temps: &Vec<Float>, rng: &mut impl Rng
    ) -> (Vec<TimeFloat>, Vec<Float>) 
    {

            let t_rand: TimeFloat = generate_open_random(0.0,1.0,rng) * self.time_birth;
            let temp_rand:Float = generate_open_random(-0.5, 0.5, rng) * self.temp_birth;

            let new_time = current_times[k] - t_rand * (current_times[k]-current_times[k+1]).abs();
            let new_temp = current_temps[k] - t_rand * (current_temps[k]-current_temps[k+1]) - temp_rand;
            
            let mut times = current_times.clone();
            let mut temps = current_temps.clone();
            times.insert(k+1, new_time);
            temps.insert(k+1, new_temp);
            
            (times, temps)
    }
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
pub struct RandomWalkConstructor {
    pub age_range: (TimeFloat, TimeFloat),
    pub temp_range: (Float, Float),
    pub fixed_point_last: bool,
    pub maintain_constraint_boxes: bool,
    pub monotonicity: char,
    pub move_probs: MoveProbabilities,
    pub move_sigmas: MoveSigmas,
}

impl RandomWalkConstructor{

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
    pub fn update_self(&mut self, new: RandomWalkConstructor) {
        self.age_range = new.age_range;
        self.temp_range= new.temp_range;
        self.fixed_point_last = new.fixed_point_last;
        self.maintain_constraint_boxes = new.maintain_constraint_boxes;
        self.monotonicity = new.monotonicity;
        self.move_probs = new.move_probs;
        self.move_sigmas = new.move_sigmas;
    }

    fn pick_k(&self, k:usize, rng: &mut impl Rng) -> usize
    {
        if self.fixed_point_last{
            generate_closed_random(0, k, rng)
        } else {
            generate_closed_random(1, k+1, rng)
        }
    }

    pub fn construct_new_profile(&self, current_profile: &TimeTempProfile,  rng: &mut impl Rng
    )-> Result<(TimeTempProfile,MoveKind,usize), String> 
    {
        let move_type = self.move_probs.choose_step(rng, current_profile.k());
        let k_choice = match move_type {
            MoveKind::Birth => {
                let time_gaps: Vec<Float> = current_profile.times
                                                    .windows(2)
                                                    .map(|pair| (pair[1] - pair[0]) .abs())
                                                    .collect();

                generate_weighted_index_random(&time_gaps,rng)
            },
            _ => self.pick_k(current_profile.k(), rng)
        };

        let (times, temps) = match move_type {

            MoveKind::Birth => {
                self.move_sigmas.birth(k_choice, &current_profile.times, &current_profile.temps, rng)
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
        Ok((TimeTempProfile{times,temps},move_type, k_choice))

    }
}


pub enum ProfileType{
    Random {constructor: RandomProfileConstructor},
    Walk   {constructor: RandomProfileConstructor, walker: RandomWalkConstructor}
}

impl ProfileType {
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
