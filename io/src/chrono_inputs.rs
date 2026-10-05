// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only
use crate::errors::InputError;
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use std::fs::read_to_string;
use std::path::Path;

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct ChronSetup {
    pub fill_ratio: Vec<Float>,
    pub fill_ratio_error: Vec<Float>,
    pub random_paths: bool,
}

impl Default for ChronSetup {
    fn default() -> Self {
        Self {
            fill_ratio: vec![0.0],
            fill_ratio_error: vec![0.1],
            random_paths: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct RandomWalkInputs {
    pub time_unit: TimeUnit,
    pub temp_unit: TemperatureUnit,
    pub step_max: usize,
    pub time_start_range: Vec<TimeFloat>,
    pub time_end_range: Vec<TimeFloat>,
    pub temp_start_range: Vec<Float>,
    pub temp_end_range: Vec<Float>,
    pub temp_range: Vec<Float>,
    pub monotonic: bool,
}

impl Default for RandomWalkInputs {
    fn default() -> Self {
        Self {
            time_unit: TimeUnit::Second,
            temp_unit: TemperatureUnit::Celsius,
            step_max: 100,
            time_start_range: vec![0.0, 0.0],
            time_end_range: vec![160.0, 160.0],
            temp_start_range: vec![0.0, 0.0],
            temp_end_range: vec![800.0, 800.0],
            temp_range: vec![0.0, 800.0],
            monotonic: true,
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
    pub random_walks: RandomWalkInputs,
    pub grid_options: GridOptions,
}
impl Default for ChronologyInputs {
    fn default() -> Self {
        Self {
            setup: ChronSetup::default(),
            random_walks: RandomWalkInputs::default(),
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
