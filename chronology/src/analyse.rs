// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only
//! Convert stored time-temperature profiles into a two-dimensional occupancy grid.
use crate::profiles::TimeTempProfile;

use common::numeric::{Float, Numeric, TimeFloat};
use io::chrono_inputs::{GridOptions};
use io::chrono_outputs::{read_all_chrono_batches, ProfileGrid};
use io::chrono_plotting::{
    plot_profile_grid_heatmap as plot_grid_heatmap, ProfileHeatmapOptions, ProfilePath,
};
use io::errors::PlotError;

use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds<T: Numeric> {
    pub min: T,
    pub max: T,
}

impl<T: Numeric> Bounds<T> {
    fn new(a: T, b: T, name: &str) -> Result<Self, String> {
        let min = a.min(b);
        let max = a.max(b);

        if !min.is_finite() || !max.is_finite() || min >= max {
            return Err(format!(
                "{name} bounds must be finite and have non-zero width"
            ));
        }

        Ok(Self { min, max })
    }

    fn include(&mut self, value: T) {
        self.min = self.min.min(value);
        self.max = self.max.max(value);
    }
}

fn validate_profile(profile: &TimeTempProfile) -> Result<(), String> {
    if profile.times.len() != profile.temps.len() {
        return Err(format!(
            "profile has {} time values but {} temperature values",
            profile.times.len(),
            profile.temps.len(),
        ));
    }

    if profile.times.len() < 2 {
        return Err("profile must contain at least two time-temperature points".to_string());
    }

    if profile.times.iter().any(|value| !value.is_finite())
        || profile.temps.iter().any(|value| !value.is_finite())
    {
        return Err("profile contains a non-finite value".to_string());
    }

    Ok(())
}

fn update_detected_bounds<T: Numeric>(values: &[T], bounds: &mut Option<Bounds<T>>) {
    for &value in values {
        match bounds {
            Some(bounds) => bounds.include(value),
            None => {
                *bounds = Some(Bounds {
                    min: value,
                    max: value,
                });
            }
        }
    }
}

/// Return the supplied bounds, or scan the stored batches for missing bounds.
pub fn find_bounds(
    input_path: impl AsRef<Path>,
    time_bounds: Option<(TimeFloat, TimeFloat)>,
    temperature_bounds: Option<(Float, Float)>,
) -> Result<(Bounds<TimeFloat>, Bounds<Float>), String> {
    let find_time = time_bounds.is_none();
    let find_temperature = temperature_bounds.is_none();
    let mut detected_time = None;
    let mut detected_temperature = None;

    if find_time || find_temperature {
        let reader = read_all_chrono_batches::<TimeTempProfile>(input_path)
            .map_err(|error| error.to_string())?;
        for batch in reader {
            for profile in batch.map_err(|error| error.to_string())? {
                validate_profile(&profile)?;
                if find_time {
                    update_detected_bounds(&profile.times, &mut detected_time);
                }
                if find_temperature {
                    update_detected_bounds(&profile.temps, &mut detected_temperature);
                }
            }
        }
    }

    let time_bounds = match time_bounds {
        Some((start, end)) => Bounds::new(start, end, "time")?,
        None => {
            let bounds = detected_time
                .ok_or_else(|| "no time values were found in the profile data".to_string())?;
            Bounds::new(bounds.min, bounds.max, "time")?
        }
    };

    let temperature_bounds = match temperature_bounds {
        Some((start, end)) => Bounds::new(start, end, "temperature")?,
        None => {
            let bounds = detected_temperature.ok_or_else(|| {
                "no temperature values were found in the profile data".to_string()
            })?;
            Bounds::new(bounds.min, bounds.max, "temperature")?
        }
    };

    Ok((time_bounds, temperature_bounds))
}

/// Clip a line segment to the grid rectangle using the Liang-Barsky algorithm.
fn clip_segment(
    time0: TimeFloat,
    temperature0: Float,
    time1: TimeFloat,
    temperature1: Float,
    time_bounds: Bounds<TimeFloat>,
    temperature_bounds: Bounds<Float>,
) -> Option<(TimeFloat, Float, TimeFloat, Float)> {
    let delta_time = time1 - time0;
    let delta_temperature = temperature1 - temperature0;
    let mut enter: TimeFloat = 0.0;
    let mut leave: TimeFloat = 1.0;

    let tests = [
        (-delta_time, time0 - time_bounds.min),
        (delta_time, time_bounds.max - time0),
        (-delta_temperature, temperature0 - temperature_bounds.min),
        (delta_temperature, temperature_bounds.max - temperature0),
    ];

    for (direction, distance) in tests {
        if direction == 0.0 {
            if distance < 0.0 {
                return None;
            }
            continue;
        }

        let parameter = distance / direction;
        if direction < 0.0 {
            enter = enter.max(parameter);
        } else {
            leave = leave.min(parameter);
        }

        if enter > leave {
            return None;
        }
    }

    Some((
        time0 + enter * delta_time,
        temperature0 + enter * delta_temperature,
        time0 + leave * delta_time,
        temperature0 + leave * delta_temperature,
    ))
}

fn cell_index<T: Numeric>(value: T, bounds: Bounds<T>, cells: usize) -> Option<usize> {
    let value = value.to_float();
    let min = bounds.min.to_float();
    let max = bounds.max.to_float();

    if value < min || value > max {
        return None;
    }

    let index = (((value - min) / (max - min)) * cells as Float).floor() as usize;
    Some(index.min(cells - 1))
}

fn mark_point(
    time: TimeFloat,
    temperature: Float,
    time_bounds: Bounds<TimeFloat>,
    temperature_bounds: Bounds<Float>,
    time_cells: usize,
    temperature_cells: usize,
    profile_number: u64,
    visited_by_profile: &mut [u64],
    weights: &mut [u64],
) {
    let Some(time_index) = cell_index(time, time_bounds, time_cells) else {
        return;
    };
    let Some(temperature_index) = cell_index(temperature, temperature_bounds, temperature_cells)
    else {
        return;
    };

    let index = temperature_index * time_cells + time_index;
    if visited_by_profile[index] != profile_number {
        visited_by_profile[index] = profile_number;
        weights[index] += 1;
    }
}

fn same_parameter(left: TimeFloat, right: TimeFloat) -> bool {
    let scale = left.abs().max(right.abs()).max(1.0);
    (left - right).abs() <= 32.0 * TimeFloat::EPSILON * scale
}

/// Mark every cell whose interior is crossed by one profile segment.
///
/// Cells are half-open on their upper edges, except that the overall maximum
/// belongs to the final cell. A segment lying exactly on an internal boundary
/// is assigned to the cell above/right of that boundary.

fn rasterize_segment(
    time0: TimeFloat,
    temperature0: Float,
    time1: TimeFloat,
    temperature1: Float,
    time_bounds: Bounds<TimeFloat>,
    temperature_bounds: Bounds<Float>,
    time_cells: usize,
    temperature_cells: usize,
    profile_number: u64,
    visited_by_profile: &mut [u64],
    weights: &mut [u64],
) {
    let Some((time0, temperature0, time1, temperature1)) = clip_segment(
        time0,
        temperature0,
        time1,
        temperature1,
        time_bounds,
        temperature_bounds,
    ) else {
        return;
    };

    let delta_time = time1 - time0;
    let delta_temperature = temperature1 - temperature0;

    if delta_time == 0.0 && delta_temperature == 0.0 {
        mark_point(
            time0,
            temperature0,
            time_bounds,
            temperature_bounds,
            time_cells,
            temperature_cells,
            profile_number,
            visited_by_profile,
            weights,
        );
        return;
    }

    let mut crossings = Vec::with_capacity(time_cells + temperature_cells + 2);
    crossings.push(0.0);
    crossings.push(1.0);

    if delta_time != 0.0 {
        let width = (time_bounds.max - time_bounds.min) / time_cells as TimeFloat;
        for boundary_index in 1..time_cells {
            let boundary = time_bounds.min + boundary_index as TimeFloat * width;
            let parameter = (boundary - time0) / delta_time;
            if parameter > 0.0 && parameter < 1.0 {
                crossings.push(parameter);
            }
        }
    }

    if delta_temperature != 0.0 {
        let height = (temperature_bounds.max - temperature_bounds.min) / temperature_cells as Float;
        for boundary_index in 1..temperature_cells {
            let boundary = temperature_bounds.min + boundary_index as Float * height;
            let parameter = (boundary - temperature0) / delta_temperature;
            if parameter > 0.0 && parameter < 1.0 {
                crossings.push(parameter);
            }
        }
    }

    crossings.sort_by(TimeFloat::total_cmp);
    crossings.dedup_by(|left, right| same_parameter(*left, *right));

    for interval in crossings.windows(2) {
        let midpoint = (interval[0] + interval[1]) * 0.5;
        mark_point(
            time0 + midpoint * delta_time,
            temperature0 + midpoint * delta_temperature,
            time_bounds,
            temperature_bounds,
            time_cells,
            temperature_cells,
            profile_number,
            visited_by_profile,
            weights,
        );
    }
}

fn add_profile(
    profile: &TimeTempProfile,
    grid: &mut ProfileGrid,
    visited_by_profile: &mut [u64],
) -> Result<(), String> {
    validate_profile(profile)?;

    let times_are_ordered = if grid.time_unit.is_ka_or_ma() {
        profile.times.windows(2).all(|times| times[0] > times[1])
    } else {
        profile.times.windows(2).all(|times| times[0] < times[1])
    };

    if !times_are_ordered {
        let direction = if grid.time_unit.is_ka_or_ma() {
            "strictly decreasing"
        } else {
            "strictly increasing"
        };
        return Err(format!(
            "profile times must be {direction} for {}",
            grid.time_unit
        ));
    }

    let time_bounds = Bounds::new(grid.time_bounds.0, grid.time_bounds.1, "time")?;
    let temperature_bounds = Bounds::new(
        grid.temperature_bounds.0,
        grid.temperature_bounds.1,
        "temperature",
    )?;

    for (times, temperatures) in profile.times.windows(2).zip(profile.temps.windows(2)) {
        rasterize_segment(
            times[0],
            temperatures[0],
            times[1],
            temperatures[1],
            time_bounds,
            temperature_bounds,
            grid.time_cells,
            grid.temperature_cells,
            grid.accepted_profiles,
            visited_by_profile,
            &mut grid.weights,
        );
    }

    Ok(())
}

/// Read all stored profiles and accumulate their paths into an evenly spaced grid.
pub fn build_profile_grid(
    input_path: impl AsRef<Path>,
    options: &GridOptions,
) -> Result<ProfileGrid, String> {
    if options.time_cells == 0 || options.temperature_cells == 0 {
        return Err("grid dimensions must be greater than zero".to_string());
    }

    let (time_bounds, temperature_bounds) =
        find_bounds(&input_path, options.time_bounds, options.temperature_bounds)?;

    let reader = read_all_chrono_batches::<TimeTempProfile>(&input_path)
        .map_err(|error| error.to_string())?;
    let time_unit = reader.time_unit;
    let temperature_unit = reader.temp_unit;

    let cell_count = options
        .time_cells
        .checked_mul(options.temperature_cells)
        .ok_or_else(|| "grid dimensions are too large".to_string())?;

    let mut grid = ProfileGrid {
        accepted_profiles: 0,
        time_bounds: (time_bounds.min, time_bounds.max),
        temperature_bounds: (temperature_bounds.min, temperature_bounds.max),
        time_unit,
        temperature_unit,
        time_cells: options.time_cells,
        temperature_cells: options.temperature_cells,
        weights: vec![0; cell_count],
    };

    let mut visited_by_profile = vec![0_u64; cell_count];

    for batch in reader {
        for profile in batch.map_err(|error| error.to_string())? {
            grid.accepted_profiles = grid
                .accepted_profiles
                .checked_add(1)
                .ok_or_else(|| "accepted profile count overflowed u64".to_string())?;
            add_profile(&profile, &mut grid, &mut visited_by_profile)?;
        }
    }

    Ok(grid)
}

/// Plot a profile-grid CSV with an optional known time-temperature profile.
///
/// The known profile must use the same units as the CSV.
pub fn plot_profile_grid_heatmap(
    csv_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    options: &ProfileHeatmapOptions,
    real_profile: Option<&TimeTempProfile>,
) -> Result<(), PlotError> {
    let overlay = real_profile.map(|profile| ProfilePath {
        times: &profile.times,
        temperatures: &profile.temps,
    });
    plot_grid_heatmap(csv_path, output_path, options, overlay)
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::constants::temperature::TemperatureUnit;
    use common::constants::time::TimeUnit;
    use io::chrono_outputs::{
        append_chronology_profile_to_experiment_file, write_temporary_chronology_units,
    };
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_output_path(label: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_analysis_{label}_{}_{}.bin.gz",
            std::process::id(),
            unique,
        ))
    }

    fn write_profiles(path: &Path, profiles: &[TimeTempProfile]) {
        write_temporary_chronology_units(path, TimeUnit::Second, TemperatureUnit::Celsius).unwrap();
        append_chronology_profile_to_experiment_file(path, profiles).unwrap();
    }

    #[test]
    fn detects_only_bounds_that_were_not_supplied() {
        let path = temporary_output_path("bounds");
        write_profiles(
            &path,
            &[TimeTempProfile {
                times: vec![0.0, 5.0, 10.0],
                temps: vec![-2.0, 4.0, 8.0],
            }],
        );

        let (time, temperature) = find_bounds(&path, Some((20.0, 0.0)), None).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(
            time,
            Bounds {
                min: 0.0,
                max: 20.0
            }
        );
        assert_eq!(
            temperature,
            Bounds {
                min: -2.0,
                max: 8.0
            }
        );
    }

    #[test]
    fn diagonal_marks_only_the_two_cells_it_crosses() {
        let path = temporary_output_path("diagonal");
        write_profiles(
            &path,
            &[TimeTempProfile {
                times: vec![0.0, 2.0],
                temps: vec![0.0, 2.0],
            }],
        );

        let grid = build_profile_grid(
            &path,
            &GridOptions {
                time_cells: 2,
                temperature_cells: 2,
                time_bounds: None,
                temperature_bounds: None,
            },
        )
        .unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(grid.accepted_profiles, 1);
        assert_eq!(grid.weights, vec![1, 0, 0, 1]);
    }

    #[test]
    fn a_profile_contributes_at_most_once_to_each_cell() {
        let path = temporary_output_path("deduplicate");
        write_profiles(
            &path,
            &[
                TimeTempProfile {
                    times: vec![0.0, 1.0, 2.0],
                    temps: vec![0.25, 0.75, 0.25],
                },
                TimeTempProfile {
                    times: vec![0.0, 2.0],
                    temps: vec![0.5, 0.5],
                },
            ],
        );

        let grid = build_profile_grid(
            &path,
            &GridOptions {
                time_cells: 1,
                temperature_cells: 1,
                time_bounds: Some((0.0, 2.0)),
                temperature_bounds: Some((0.0, 1.0)),
            },
        )
        .unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(grid.accepted_profiles, 2);
        assert_eq!(grid.weights, vec![2]);
    }

    #[test]
    fn specified_bounds_clip_segments_without_out_of_range_indices() {
        let path = temporary_output_path("clip");
        write_profiles(
            &path,
            &[TimeTempProfile {
                times: vec![-1.0, 3.0],
                temps: vec![-1.0, 3.0],
            }],
        );

        let grid = build_profile_grid(
            &path,
            &GridOptions {
                time_cells: 2,
                temperature_cells: 2,
                time_bounds: Some((0.0, 2.0)),
                temperature_bounds: Some((0.0, 2.0)),
            },
        )
        .unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(grid.weights, vec![1, 0, 0, 1]);
    }

    #[test]
    fn rejects_time_order_that_disagrees_with_the_unit() {
        let path = temporary_output_path("direction");
        write_profiles(
            &path,
            &[TimeTempProfile {
                times: vec![2.0, 0.0],
                temps: vec![0.0, 1.0],
            }],
        );

        let error = build_profile_grid(
            &path,
            &GridOptions {
                time_cells: 2,
                temperature_cells: 2,
                time_bounds: Some((0.0, 2.0)),
                temperature_bounds: Some((0.0, 1.0)),
            },
        )
        .err()
        .expect("decreasing seconds should be rejected");
        fs::remove_file(path).unwrap();

        assert!(error.contains("strictly increasing"));
    }

    #[test]
    fn geological_time_accepts_decreasing_profiles() {
        let path = temporary_output_path("geological");
        write_temporary_chronology_units(&path, TimeUnit::KAnnum, TemperatureUnit::Celsius)
            .unwrap();
        append_chronology_profile_to_experiment_file(
            &path,
            &[TimeTempProfile {
                times: vec![2.0, 0.0],
                temps: vec![0.0, 2.0],
            }],
        )
        .unwrap();

        let grid = build_profile_grid(
            &path,
            &GridOptions {
                time_cells: 2,
                temperature_cells: 2,
                time_bounds: None,
                temperature_bounds: None,
            },
        )
        .unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(grid.time_bounds, (0.0, 2.0));
        assert_eq!(grid.weights, vec![0, 1, 1, 0]);
    }
}
