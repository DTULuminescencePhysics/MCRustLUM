// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Load a chronology profile-grid CSV and plot its weighted path heatmap.

use crate::errors::PlotError;
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use plotters::prelude::*;
use std::path::Path;
use std::str::FromStr;

/// One lower/upper weighted-quantile interval drawn around the median path.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantileBand {
    pub lower: Float,
    pub upper: Float,
    pub label: String,
}

impl QuantileBand {
    /// Construct an arbitrary interval such as `0.10..0.90`.
    pub fn new(lower: Float, upper: Float, label: impl Into<String>) -> Result<Self, String> {
        if !lower.is_finite() || !upper.is_finite() || lower < 0.0 || upper > 1.0 || lower >= upper
        {
            return Err("quantile bounds must satisfy 0 <= lower < upper <= 1".to_string());
        }

        Ok(Self {
            lower,
            upper,
            label: label.into(),
        })
    }

    /// Construct a central coverage interval.
    ///
    /// For example, `central(0.90)` produces the 5th-to-95th percentile band,
    /// while `central(0.60)` produces the 20th-to-80th percentile band.
    pub fn central(coverage: Float) -> Result<Self, String> {
        if !coverage.is_finite() || coverage <= 0.0 || coverage >= 1.0 {
            return Err("central coverage must be finite and between 0 and 1".to_string());
        }

        let tail = (1.0 - coverage) * 0.5;
        Self::new(
            tail,
            1.0 - tail,
            format!("{:.0}% interval", coverage * 100.0),
        )
    }
}

/// Rendering controls for a chronology heatmap.
#[derive(Debug, Clone)]
pub struct ProfileHeatmapOptions {
    pub width: u32,
    pub height: u32,
    pub caption: String,
    /// Intervals are drawn in this order, so place wider bands first.
    pub bands: Vec<QuantileBand>,
    /// Use `ln(1 + weight)` to make low-density paths more visible.
    pub logarithmic_colours: bool,
}

impl Default for ProfileHeatmapOptions {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 800,
            caption: "Time-temperature path heatmap".to_string(),
            bands: vec![
                QuantileBand::central(0.90).expect("constant coverage is valid"),
                QuantileBand::central(0.60).expect("constant coverage is valid"),
            ],
            logarithmic_colours: false,
        }
    }
}

/// Optional known time-temperature path overlaid on the inferred paths.
#[derive(Debug, Clone, Copy)]
pub struct ProfilePath<'a> {
    pub times: &'a [TimeFloat],
    pub temperatures: &'a [Float],
}

/// Parsed contents of `write_profile_grid_csv` output.
#[derive(Debug, Clone)]
pub struct ProfileGridCsv {
    pub accepted_profiles: u64,
    /// Start and end are retained in display order (`0 -> t` or `t -> 0`).
    pub time_start: TimeFloat,
    pub time_end: TimeFloat,
    pub time_unit: TimeUnit,
    pub temperature_min: Float,
    pub temperature_max: Float,
    pub temperature_unit: TemperatureUnit,
    /// Rows are stored high-temperature first, exactly as they appear in CSV.
    pub weights: Vec<Vec<u64>>,
}

impl ProfileGridCsv {
    pub fn time_cells(&self) -> usize {
        self.weights.first().map_or(0, Vec::len)
    }

    pub fn temperature_cells(&self) -> usize {
        self.weights.len()
    }

    fn time_width(&self) -> TimeFloat {
        (self.time_end - self.time_start) / self.time_cells() as TimeFloat
    }

    fn temperature_height(&self) -> Float {
        (self.temperature_max - self.temperature_min) / self.temperature_cells() as Float
    }

    fn time_centre(&self, column: usize) -> TimeFloat {
        self.time_start + (column as TimeFloat + 0.5) * self.time_width()
    }

    fn temperature_centre(&self, row: usize) -> Float {
        self.temperature_max - (row as Float + 0.5) * self.temperature_height()
    }

    /// Weighted temperature quantile for one displayed time column.
    fn column_quantile(&self, column: usize, quantile: Float) -> Option<Float> {
        let total = self.weights.iter().map(|row| row[column]).sum::<u64>();
        if total == 0 {
            return None;
        }

        // Iterate from low to high temperature. CSV rows have the opposite order.
        let target = quantile * total.saturating_sub(1) as Float;
        let mut cumulative = 0_u64;
        for row in (0..self.temperature_cells()).rev() {
            cumulative += self.weights[row][column];
            if cumulative as Float > target {
                return Some(self.temperature_centre(row));
            }
        }

        Some(self.temperature_centre(0))
    }
}

fn invalid(path: &Path, message: impl Into<String>) -> PlotError {
    PlotError::InvalidData {
        path: path.to_path_buf(),
        message: message.into(),
    }
}

fn parse_field<T: std::str::FromStr>(
    path: &Path,
    value: Option<&str>,
    name: &str,
) -> Result<T, PlotError> {
    value
        .ok_or_else(|| invalid(path, format!("missing {name}")))?
        .parse()
        .map_err(|_| invalid(path, format!("invalid {name}")))
}

fn parse_time_unit(path: &Path, value: Option<&str>) -> Result<TimeUnit, PlotError> {
    let value = value.ok_or_else(|| invalid(path, "missing time unit"))?;
    TimeUnit::from_str(value).map_err(|e| invalid(path, e))     
}

fn parse_temperature_unit(path: &Path, value: Option<&str>) -> Result<TemperatureUnit, PlotError> {
    let value = value.ok_or_else(|| invalid(path, "missing temperature unit"))?;
    TemperatureUnit::from_str(value).map_err(|e| invalid(path, e))     
}

/// Load and validate a CSV created by `write_profile_grid_csv`.
pub fn read_profile_grid_csv(path: impl AsRef<Path>) -> Result<ProfileGridCsv, PlotError> {
    let path = path.as_ref();
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(path)
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?;
    let mut records = reader.records();

    let count = records
        .next()
        .ok_or_else(|| invalid(path, "missing accepted-profile count"))?
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?;
    if count.len() != 1 {
        return Err(invalid(
            path,
            "the first row must contain only the accepted-profile count",
        ));
    }
    let accepted_profiles = parse_field(path, count.get(0), "accepted-profile count")?;

    let time = records
        .next()
        .ok_or_else(|| invalid(path, "missing time metadata row"))?
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?;
    if time.len() != 3 {
        return Err(invalid(
            path,
            "the time row must contain start, end, and unit",
        ));
    }
    let time_start: TimeFloat = parse_field(path, time.get(0), "start time")?;
    let time_end: TimeFloat = parse_field(path, time.get(1), "end time")?;
    let time_unit  =  parse_time_unit(path, time.get(2))?;
    
    
    let temperature = records
        .next()
        .ok_or_else(|| invalid(path, "missing temperature metadata row"))?
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?;
    if temperature.len() != 3 {
        return Err(invalid(
            path,
            "the temperature row must contain minimum, maximum, and unit",
        ));
    }
    let temperature_min: Float = parse_field(path, temperature.get(0), "minimum temperature")?;
    let temperature_max: Float = parse_field(path, temperature.get(1), "maximum temperature")?;
    let temperature_unit = parse_temperature_unit(path, temperature.get(2))?;

    if !time_start.is_finite() || !time_end.is_finite() || time_start == time_end {
        return Err(invalid(path, "time bounds must be finite and distinct"));
    }
    if !temperature_min.is_finite()
        || !temperature_max.is_finite()
        || temperature_min >= temperature_max
    {
        return Err(invalid(
            path,
            "temperature bounds must be finite and increasing",
        ));
    }
    let expected_reverse = time_unit.is_ka_or_ma();
    if (time_start > time_end) != expected_reverse {
        return Err(invalid(
            path,
            "time direction does not agree with the stored time unit",
        ));
    }

    let mut weights = Vec::new();
    let mut columns = None;
    for record in records {
        let record = record.map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?;
        if record.is_empty() {
            continue;
        }
        match columns {
            Some(columns) if record.len() != columns => {
                return Err(invalid(path, "grid rows have different lengths"));
            }
            None => columns = Some(record.len()),
            _ => {}
        }
        let row = record
            .iter()
            .map(|value| {
                value
                    .parse::<u64>()
                    .map_err(|_| invalid(path, format!("invalid grid weight {value:?}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        weights.push(row);
    }

    if weights.is_empty() || columns == Some(0) {
        return Err(invalid(path, "profile grid is empty"));
    }

    Ok(ProfileGridCsv {
        accepted_profiles,
        time_start,
        time_end,
        time_unit,
        temperature_min,
        temperature_max,
        temperature_unit,
        weights,
    })
}


fn draw_error(path: &Path, error: impl std::fmt::Debug) -> PlotError {
    PlotError::Draw {
        path: path.to_path_buf(),
        message: format!("{error:?}"),
    }
}

fn validate_options(output: &Path, options: &ProfileHeatmapOptions) -> Result<(), PlotError> {
    if options.width == 0 || options.height == 0 {
        return Err(invalid(output, "plot width and height must be non-zero"));
    }
    for band in &options.bands {
        QuantileBand::new(band.lower, band.upper, band.label.clone())
            .map_err(|message| invalid(output, message))?;
    }
    Ok(())
}

fn validate_overlay(
    output: &Path,
    path: ProfilePath<'_>,
    time_unit: TimeUnit,
) -> Result<(), PlotError> {
    if path.times.len() != path.temperatures.len() || path.times.len() < 2 {
        return Err(invalid(
            output,
            "real profile must contain matching time and temperature vectors with at least two points",
        ));
    }
    if path.times.iter().any(|value| !value.is_finite())
        || path.temperatures.iter().any(|value| !value.is_finite())
    {
        return Err(invalid(output, "real profile contains a non-finite value"));
    }
    let correctly_ordered = if time_unit.is_ka_or_ma() {
        path.times.windows(2).all(|times| times[0] > times[1])
    } else {
        path.times.windows(2).all(|times| times[0] < times[1])
    };
    if !correctly_ordered {
        return Err(invalid(
            output,
            "real-profile time direction does not agree with the grid time unit",
        ));
    }
    Ok(())
}

/// Plot a heatmap, median path, quantile bands, and optional known path.
///
/// The optional path must use the same units recorded in the grid CSV.
pub fn plot_profile_grid_heatmap(
    csv_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    options: &ProfileHeatmapOptions,
    real_profile: Option<ProfilePath<'_>>,
    
) -> Result<(), PlotError> {
    let csv_path = csv_path.as_ref();
    let output_path = output_path.as_ref();
    validate_options(output_path, options)?;
    let grid = read_profile_grid_csv(csv_path)?;
    if let Some(profile) = real_profile {
        validate_overlay(output_path, profile, grid.time_unit)?;
    }

    let root = BitMapBackend::new(output_path, (options.width, options.height)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|error| draw_error(output_path, error))?;

    let caption = format!(
        "{} ({} accepted profiles)",
        options.caption, grid.accepted_profiles
    );
    let mut chart = ChartBuilder::on(&root)
        .caption(caption, ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(55)
        .y_label_area_size(75)
        .build_cartesian_2d(
            grid.time_start..grid.time_end,
            grid.temperature_min..grid.temperature_max,
        )
        .map_err(|error| draw_error(output_path, error))?;

    chart
        .configure_mesh()
        .x_desc(format!("Time ({})", grid.time_unit))
        .y_desc(format!(
            "Temperature ({})",
            grid.temperature_unit
        ))
        .draw()
        .map_err(|error| draw_error(output_path, error))?;

    let maximum_weight = grid.weights.iter().flatten().copied().max().unwrap_or(0);
    let time_width = grid.time_width();
    let temperature_height = grid.temperature_height();
    if maximum_weight > 0 {
        chart
            .draw_series(grid.weights.iter().enumerate().flat_map(|(row, weights)| {
                weights
                    .iter()
                    .enumerate()
                    .filter_map(move |(column, &weight)| {
                        if weight == 0 {
                            return None;
                        }
                        let intensity = if options.logarithmic_colours {
                            (1.0 + weight as Float).ln() / (1.0 + maximum_weight as Float).ln()
                        } else {
                            weight as Float / maximum_weight as Float
                        };
                        let time0 = grid.time_start + column as TimeFloat * time_width;
                        let time1 = time0 + time_width;
                        let temperature1 = grid.temperature_max - row as Float * temperature_height;
                        let temperature0 = temperature1 - temperature_height;
                        let colour = HSLColor((1.0 - intensity) * 240.0 / 360.0, 0.95, 0.52);
                        Some(Rectangle::new(
                            [(time0, temperature0), (time1, temperature1)],
                            colour.mix(0.78).filled(),
                        ))
                    })
            }))
            .map_err(|error| draw_error(output_path, error))?;
    }

    let band_colours = [
        RGBColor(40, 120, 220),
        RGBColor(0, 150, 110),
        RGBColor(160, 80, 180),
        RGBColor(220, 130, 20),
    ];
    for (band_index, band) in options.bands.iter().enumerate() {
        let colour = band_colours[band_index % band_colours.len()];
        let mut run = Vec::<(TimeFloat, Float, Float)>::new();
        let mut runs = Vec::<Vec<(TimeFloat, Float, Float)>>::new();
        for column in 0..grid.time_cells() {
            match (
                grid.column_quantile(column, band.lower),
                grid.column_quantile(column, band.upper),
            ) {
                (Some(lower), Some(upper)) => {
                    run.push((grid.time_centre(column), lower, upper));
                }
                _ if !run.is_empty() => runs.push(std::mem::take(&mut run)),
                _ => {}
            }
        }
        if !run.is_empty() {
            runs.push(run);
        }

        for (run_index, run) in runs.iter().enumerate() {
            if run.len() < 2 {
                continue;
            }
            let mut polygon = run
                .iter()
                .map(|&(time, _, upper)| (time, upper))
                .collect::<Vec<_>>();
            polygon.extend(run.iter().rev().map(|&(time, lower, _)| (time, lower)));

            let series = chart
                .draw_series(std::iter::once(Polygon::new(
                    polygon,
                    colour.mix(0.16).filled(),
                )))
                .map_err(|error| draw_error(output_path, error))?;
            if run_index == 0 {
                series.label(band.label.clone()).legend(move |(x, y)| {
                    Rectangle::new([(x, y - 5), (x + 20, y + 5)], colour.mix(0.24).filled())
                });
            }

            chart
                .draw_series(LineSeries::new(
                    run.iter().map(|&(time, lower, _)| (time, lower)),
                    colour.stroke_width(1),
                ))
                .map_err(|error| draw_error(output_path, error))?;
            chart
                .draw_series(LineSeries::new(
                    run.iter().map(|&(time, _, upper)| (time, upper)),
                    colour.stroke_width(1),
                ))
                .map_err(|error| draw_error(output_path, error))?;
        }
    }

    let mut median_run = Vec::new();
    let mut median_runs = Vec::new();
    for column in 0..grid.time_cells() {
        if let Some(median) = grid.column_quantile(column, 0.5) {
            median_run.push((grid.time_centre(column), median));
        } else if !median_run.is_empty() {
            median_runs.push(std::mem::take(&mut median_run));
        }
    }
    if !median_run.is_empty() {
        median_runs.push(median_run);
    }
    for (run_index, run) in median_runs.iter().enumerate() {
        if run.is_empty() {
            continue;
        }
        let series = chart
            .draw_series(LineSeries::new(run.iter().copied(), BLACK.stroke_width(3)))
            .map_err(|error| draw_error(output_path, error))?;
        if run_index == 0 {
            series
                .label("Median path")
                .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], BLACK.stroke_width(3)));
        }
    }

    if let Some(profile) = real_profile {
        chart
            .draw_series(LineSeries::new(
                profile
                    .times
                    .iter()
                    .copied()
                    .zip(profile.temperatures.iter().copied()),
                RED.stroke_width(3),
            ))
            .map_err(|error| draw_error(output_path, error))?
            .label("Known profile")
            .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], RED.stroke_width(3)));
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()
        .map_err(|error| draw_error(output_path, error))?;
    root.present()
        .map_err(|error| draw_error(output_path, error))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_path(extension: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_heatmap_{}_{}.{}",
            std::process::id(),
            unique,
            extension,
        ))
    }

    #[test]
    fn central_intervals_have_the_expected_quantiles() {
        let outer = QuantileBand::central(0.90).unwrap();
        let inner = QuantileBand::central(0.60).unwrap();

        assert!((outer.lower - 0.05).abs() < 1.0e-12);
        assert!((outer.upper - 0.95).abs() < 1.0e-12);
        assert!((inner.lower - 0.20).abs() < 1.0e-12);
        assert!((inner.upper - 0.80).abs() < 1.0e-12);
    }

    #[test]
    fn reads_quantiles_and_renders_reverse_time_heatmap() {
        let csv_path = temporary_path("csv");
        let png_path = temporary_path("png");
        fs::write(
            &csv_path,
            "3\n2,0,Ka annum\n0,3,Celsius\n0,1\n1,2\n2,0\n",
        )
        .unwrap();

        let grid = read_profile_grid_csv(&csv_path).unwrap();
        assert_eq!(grid.time_cells(), 2);
        assert_eq!(grid.temperature_cells(), 3);
        assert_eq!(grid.column_quantile(0, 0.5), Some(0.5));

        plot_profile_grid_heatmap(
            &csv_path,
            &png_path,
            &ProfileHeatmapOptions::default(),
            Some(ProfilePath {
                times: &[2.0, 0.0],
                temperatures: &[0.5, 2.5],
            }),
        )
        .unwrap();

        assert!(fs::metadata(&png_path).unwrap().len() > 0);
        fs::remove_file(csv_path).unwrap();
        fs::remove_file(png_path).unwrap();
    }
}
