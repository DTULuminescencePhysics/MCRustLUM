// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Chronology TOML configuration and defaults.
//!
//! Omitted fields are filled by Serde defaults. Model constructors perform
//! additional validation; loading alone does not establish a valid inference run.
use crate::errors::InputError;
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use std::fs::read_to_string;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
/// Profile-generation/acceptance modes and coordinate units for chronology inference.
///
/// Generation uses `R` (independent) or `W` (walk); acceptance uses `F`
/// (frequentist) or `B` (Bayesian). Missing TOML fields use defaults.
pub struct ChronSetup {
    /// Generation mode: `R` for independent profiles or `W` for a walk.
    pub profile_generation: char, 
    /// Statistical acceptance mode (`F` or `B`) or stored acceptance metadata.
    pub acceptance: char, 
    /// Time unit used to interpret the associated coordinates.
    pub time_unit: TimeUnit,
    /// Temperature unit associated with chronology coordinates.
    pub temp_unit: TemperatureUnit,
}

impl Default for ChronSetup {
    /// Return independent frequentist generation in seconds and Celsius.
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
/// Optional adaptive burn-in limits and target acceptance ranges for each move.
///
/// `window` is proposals per tuning window, `maximum_no_windows` caps tuning,
/// and `stable_window` requires consecutive in-range windows for all moves.
pub struct BurnIn {
    /// Whether to perform adaptive burn-in before eligible Bayesian walks.
    pub run: bool,
    /// Number of proposal attempts in one global tuning window.
    pub window: usize,
    /// Maximum number of tuning windows before production starts.
    pub maximum_no_windows: usize,
    /// Number of consecutive in-range windows required for every move kind.
    pub stable_window: usize,
    /// Inclusive target acceptance range for birth moves.
    pub birth_range:(Float, Float),
    /// Inclusive target acceptance range for death moves.
    pub death_range:(Float, Float),
    /// Inclusive target acceptance range for time-coordinate moves.
    pub time_range:(Float, Float),
    /// Inclusive target acceptance range for temperature-coordinate moves.
    pub temp_range:(Float, Float),
}

impl Default for BurnIn {
    /// Return disabled burn-in with 100-proposal windows and per-move acceptance targets.
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
/// Observed fills and ages with their measurement uncertainties.
///
/// Fill vectors follow experiment order. The acceptor currently uses only the
/// first age and age-error values. Uncertainties must be non-zero and ages
/// must use the configured profile time unit.
pub struct ChronMeasurement {
    /// Measured occupied-trap fraction for each experiment, in experiment order.
    pub fill_ratio: Vec<Float>,
    /// Non-zero measurement uncertainty for each corresponding fill ratio.
    pub fill_ratio_error: Vec<Float>,
    /// Measured profile duration in the setup time unit; only the first array entry is used.
    pub age: Vec<Float>,
    /// Non-zero age uncertainty in the setup time unit; only the first array entry is used.
    pub age_error: Vec<Float>,
}

impl Default for ChronMeasurement {
    /// Return one zero fill and age measurement, each with uncertainty 0.1.
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
/// Global profile bounds, fixed time-zero temperature range, and minimum node spacing.
///
/// Times and temperatures use the setup units. Validation requires adjacent
/// changes strictly greater than the minima and segment counts below `step_max`.
pub struct ChronConstraints {
    /// Allowed temperature range at the fixed time-zero endpoint.
    pub fixed_temp_point_range: (Float, Float),
    /// Allowed bounds for every profile time coordinate, in the setup time unit.
    pub age_range: (TimeFloat, TimeFloat),
    /// Allowed temperature bounds in profile units.
    pub temp_range: (Float, Float),
    /// Exclusive upper bound on the number of profile segments accepted by validation.
    pub step_max: usize,
    /// Strict minimum directed time difference between adjacent nodes.
    pub min_time: TimeFloat,
    /// Strict minimum absolute temperature difference between adjacent nodes.
    pub min_temp: Float,
}

impl Default for ChronConstraints {
    /// Return default global bounds and strictly positive node-spacing requirements.
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
/// Constraint boxes and subdivision rules for independent thermal histories.
///
/// Each box is `((time_min, time_max), (temperature_min, temperature_max))`.
/// Sampling codes `I`, `E`, `G` average 1, 2, or 4 uniform draws respectively.
/// Monotonicity codes are `I` (increasing), `D` (decreasing), and `V` (unrestricted).
pub struct RandomProfiles {
    /// Number of sampled constraint boxes, excluding the added fixed endpoint.
    pub constraint_box_num: usize,
    /// Constraint rectangles as paired time bounds and temperature bounds.
    pub constraint_boxes: Vec< ((TimeFloat,TimeFloat), (Float,Float))>,
    /// Recursive subdivision depth; each interval receives `2^depth - 1` interior nodes.
    pub intermediate_node_halving_num:  u32,
    /// Sampling mode: `I`, `E`, or `G` averages 1, 2, or 4 uniform draws.
    pub random_segment_kind: char, 
    /// Per-segment modes: `I` increasing, `D` decreasing, or `V` unrestricted.
    pub monotonicity: Vec<char>,
}


impl Default for RandomProfiles {
    /// Return one constraint box with decreasing-temperature subdivision settings.
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
/// Move weights, proposal widths, and constraints for the profile walk.
///
/// `p_*` values are normalized when the walker is built. `s_*` widths use
/// profile coordinate units. `k_max` limits the segment count for birth moves.
pub struct RandomWalkInputs {
    /// Maximum segment count at which further birth moves are disabled.
    pub k_max: usize,
    /// Relative selection weight for inserting a node.
    pub p_birth: Float,
    /// Relative selection weight for removing a node.
    pub p_death: Float,
    /// Relative selection weight for perturbing a time coordinate.
    pub p_time:  Float, 
    /// Relative selection weight for perturbing a temperature coordinate.
    pub p_temp:  Float, 
    /// Time interpolation scale for birth proposals.
    pub s_time_birth: TimeFloat,
    /// Temperature displacement width for birth proposals, in profile units.
    pub s_temp_birth: Float,
    /// Gaussian time-perturbation standard deviation, in profile units.
    pub s_time: TimeFloat, 
    /// Gaussian temperature-perturbation standard deviation, in profile units.
    pub s_temp: Float,
    /// Whether walk proposals must continue to satisfy box constraints.
    pub maintain_constraint_boxes: bool,
    /// Temperature ordering mode: `I` increasing, `D` decreasing, or `V` unrestricted.
    pub monotonicity: char,
}

impl Default for RandomWalkInputs {
    /// Return default birth/death/time/temperature weights and proposal scales.
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
/// Resolution and optional physical bounds of the accepted-profile occupancy grid.
///
/// Unspecified bounds are detected from stored profiles. Both cell counts
/// must be positive and explicit bounds must have finite, non-zero width.
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

    /// Return a 50 by 50 grid with bounds detected from stored profiles.
    fn default() -> Self {
        Self { time_cells: 50,
             temperature_cells: 50, 
             time_bounds: None, 
             temperature_bounds: None}
    }
}
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
/// Complete chronology configuration grouped by proposal, measurement, tuning, and grid settings.
///
/// Serde fills omitted sections and fields from their `Default` implementations.
/// Loading does not perform the runtime validation done by model constructors.
pub struct ChronologyInputs {
    /// Generation/acceptance mode selection and coordinate units.
    pub setup: ChronSetup,
    /// Observed fill and age targets and their uncertainties.
    pub measurement: ChronMeasurement,
    /// Global coordinate bounds and minimum node separation.
    pub constraints: ChronConstraints,
    /// Independent-profile constraint boxes and subdivision settings.
    pub random_profiles: RandomProfiles,
    /// Move weights and proposal widths for walk generation.
    pub random_walk: RandomWalkInputs,
    /// Optional adaptive proposal tuning settings.
    pub burn_in: BurnIn,
    /// Resolution and bounds for the final accepted-profile heatmap.
    pub grid_options: GridOptions,
}
impl Default for ChronologyInputs {
    /// Return the default value of every chronology configuration section.
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

/// Read chronology settings from TOML, filling omitted fields from defaults.
///
/// File-read and TOML-parse errors retain the source path. Physical and
/// proposal-setting validation occurs when inference state is constructed.
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
