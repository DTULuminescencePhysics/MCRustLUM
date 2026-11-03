// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only
use crate::errors::InputError;
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use std::fs::read_to_string;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct ChronSetup {
    pub profile_generation: char, 
    pub acceptance: char, 
    pub time_unit: TimeUnit,
    pub temp_unit: TemperatureUnit,
}

impl Default for ChronSetup {
    fn default() -> Self {
        Self {
            profile_generation: 'R',
            acceptance: 'F',
            time_unit: TimeUnit::Second,
            temp_unit: TemperatureUnit::Celsius,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct BurnIn {
    pub run: bool,
    pub window: usize,
    pub maximum_no_windows: usize,
    pub stable_window: usize,
    pub birth_range:(Float, Float),
    pub death_range:(Float, Float),
    pub time_range:(Float, Float),
    pub temp_range:(Float, Float),
}

impl Default for BurnIn {
    fn default() -> Self {
        Self {
            run: false,
            window: 100,
            maximum_no_windows: 1000,
            stable_window: 3,
            birth_range: (0.10, 0.40),
            death_range: (0.10, 0.40),
            time_range: (0.20, 0.50),
            temp_range: (0.20, 0.50),
        }
    }
}


#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct ChronMeasurement {
    pub fill_ratio: Vec<Float>,
    pub fill_ratio_error: Vec<Float>,
    pub age: Vec<Float>,
    pub age_error: Vec<Float>,
}

impl Default for ChronMeasurement {
    fn default() -> Self {
        Self {
            fill_ratio: vec![0.0],
            fill_ratio_error: vec![0.1],
            age: vec![0.0],
            age_error: vec![0.1],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct ChronConstraints {
    pub fixed_temp_point_range: (Float, Float),
    pub age_range: (TimeFloat, TimeFloat),
    pub temp_range: (Float, Float),
    pub step_max: usize,
    pub min_time: TimeFloat,
    pub min_temp: Float,
}

impl Default for ChronConstraints {
    fn default() -> Self {
        Self {
            fixed_temp_point_range: (0.0,5.0),
            age_range: (0.0,160.0),
            temp_range: (0.0, 160.0),
            step_max: 100,
            min_time: 1.0,
            min_temp: 1.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct RandomProfiles {
    pub constraint_box_num: usize,
    pub constraint_boxes: Vec< ((TimeFloat,TimeFloat), (Float,Float))>,
    pub intermediate_node_halving_num:  u32,
    pub random_segment_kind: char, 
    pub monotonicity: Vec<char>,
}


impl Default for RandomProfiles {
    fn default() -> Self {
        Self {
            constraint_box_num: 1,
            constraint_boxes: vec![((155.0,165.0),(155.0,165.0))],
            intermediate_node_halving_num:  10,
            random_segment_kind: 'E', 
            monotonicity: vec!['D'],
        }
    }
}


#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct RandomWalkInputs {
    pub k_max: usize,
    pub p_birth: Float,
    pub p_death: Float,
    pub p_time:  Float, 
    pub p_temp:  Float, 
    pub s_time_birth: TimeFloat,
    pub s_temp_birth: Float,
    pub s_time: TimeFloat, 
    pub s_temp: Float,
    pub maintain_constraint_boxes: bool,
    pub monotonicity: char,
}

impl Default for RandomWalkInputs {
    fn default() -> Self {
        Self {
            k_max: 100,
            p_birth: 0.15,
            p_death: 0.15,
            p_time:  0.35, 
            p_temp:  0.35, 
            s_time_birth: 20.0,
            s_temp_birth: 15.0,
            s_time: 15.0, 
            s_temp: 20.0,
            maintain_constraint_boxes: false,
            monotonicity: 'D',
        }
    }
}

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct GridOptions {
    /// Number of columns along the time axis.
    pub time_cells: usize,

    /// Number of rows along the temperature axis.
    pub temperature_cells: usize,

    /// `None` means determine the range from the stored profiles.
    pub time_bounds: Option<(TimeFloat, TimeFloat)>,

    /// `None` means determine the range from the stored profiles.
    pub temperature_bounds: Option<(Float, Float)>,
}
impl Default for GridOptions {

    fn default() -> Self {
        Self { time_cells: 50,
             temperature_cells: 50, 
             time_bounds: None, 
             temperature_bounds: None}
    }
}
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct ChronologyInputs {
    pub setup: ChronSetup,
    pub measurement: ChronMeasurement,
    pub constraints: ChronConstraints,
    pub random_profiles: RandomProfiles,
    pub random_walk: RandomWalkInputs,
    pub burn_in: BurnIn,
    pub grid_options: GridOptions,
}
impl Default for ChronologyInputs {
    fn default() -> Self {
        Self {
            setup: ChronSetup::default(),
            measurement: ChronMeasurement::default(),
            constraints: ChronConstraints::default(),
            random_profiles: RandomProfiles::default(),
            random_walk: RandomWalkInputs::default(),
            burn_in: BurnIn::default(),
            grid_options: GridOptions::default(),

        }
    }
}

pub fn read_chrono_inputs(path: impl AsRef<Path>) -> Result<ChronologyInputs, InputError> {
    let path = path.as_ref();
    let contents = read_to_string(path).map_err(|source| InputError::Read {
        path: path.to_path_buf(),
        source,
    })?;

    toml::from_str(&contents).map_err(|source| InputError::Parse {
        path: path.to_path_buf(),
        source,
    })
}
