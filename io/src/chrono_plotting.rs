// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Load a chronology profile-grid CSV and plot its weighted path heatmap.

use crate::errors::PlotError;
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use plotters::coord::Shift;
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
/// Coordinates needed to draw a stored chronology profile.
pub trait ProfileCoordinates {
    fn coordinates(&self) -> ProfilePath<'_>;
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
    pub best_fit_times: Vec<TimeFloat>,
    pub best_fit_temperatures: Vec<Float>,
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

    let best_fit_times = records
        .next()
        .ok_or_else(|| invalid(path, "missing best-fit time row"))?
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?
        .iter()
        .map(|value| {
            value
                .parse::<TimeFloat>()
                .map_err(|_| invalid(path, format!("invalid best-fit time {value:?}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let best_fit_temperatures = records
        .next()
        .ok_or_else(|| invalid(path, "missing best-fit temperature row"))?
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?
        .iter()
        .map(|value| {
            value
                .parse::<Float>()
                .map_err(|_| invalid(path, format!("invalid best-fit temperature {value:?}")))
        })
        .collect::<Result<Vec<_>, _>>()?;

    if best_fit_times.len() != best_fit_temperatures.len() || best_fit_times.len() < 2 {
        return Err(invalid(
            path,
            "best-fit path must have matching time and temperature rows with at least two values",
        ));
    }

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
    if best_fit_times.iter().any(|value| !value.is_finite())
        || best_fit_temperatures.iter().any(|value| !value.is_finite())
    {
        return Err(invalid(path, "best-fit path contains a non-finite value"));
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
        best_fit_times,
        best_fit_temperatures,
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
    validate_overlay(
        output_path,
        ProfilePath {
            times: &grid.best_fit_times,
            temperatures: &grid.best_fit_temperatures,
        },
        grid.time_unit,
    )?;
    if let Some(profile) = real_profile {
        validate_overlay(output_path, profile, grid.time_unit)?;
    }

    let root = BitMapBackend::new(output_path, (options.width, options.height)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|error| draw_error(output_path, error))?;
    let (plot_area, scale_area) = root.split_horizontally(options.width.saturating_sub(105));

    let caption = format!(
        "{} ({} accepted profiles)",
        options.caption, grid.accepted_profiles
    );
    let mut chart = ChartBuilder::on(&plot_area)
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

    let maximum_probability = 0.10;
    let time_width = grid.time_width();
    let temperature_height = grid.temperature_height();
    chart
        .draw_series(grid.weights.iter().enumerate().flat_map(|(row, weights)| {
            weights.iter().enumerate().map(move |(column, &weight)| {
                let probability = weight as Float / grid.accepted_profiles.max(1) as Float;
                let intensity = (probability / maximum_probability).clamp(0.0, 1.0);
                let time0 = grid.time_start + column as TimeFloat * time_width;
                let time1 = time0 + time_width;
                let temperature1 = grid.temperature_max - row as Float * temperature_height;
                let temperature0 = temperature1 - temperature_height;
                let colour = if weight == 0 {
                    WHITE
                } else {
                    probability_colour(intensity)
                };
                Rectangle::new(
                    [(time0, temperature0), (time1, temperature1)],
                    colour.filled(),
                )
            })
        }))
        .map_err(|error| draw_error(output_path, error))?;

    let band_colours = [
        BLACK,
        RGBColor(0, 150, 0),
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
                    colour.mix(0.06).filled(),
                )))
                .map_err(|error| draw_error(output_path, error))?;
            if run_index == 0 {
                series.label(band.label.clone()).legend(move |(x, y)| {
                    PathElement::new([(x, y), (x + 20, y)], colour.stroke_width(5))
                });
            }

            chart
                .draw_series(DashedLineSeries::new(
                    run.iter().map(|&(time, lower, _)| (time, lower)),
                    1,
                    5,
                    colour.stroke_width(5),
                ))
                .map_err(|error| draw_error(output_path, error))?;
            chart
                .draw_series(DashedLineSeries::new(
                    run.iter().map(|&(time, _, upper)| (time, upper)),
                    1,
                    5,
                    colour.stroke_width(5),
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
    let median_colour = RGBColor(255, 140, 0);
    for (run_index, run) in median_runs.iter().enumerate() {
        if run.is_empty() {
            continue;
        }
        let series = chart
            .draw_series(LineSeries::new(run.iter().copied(), median_colour.stroke_width(5)))
            .map_err(|error| draw_error(output_path, error))?;
        if run_index == 0 {
            series
                .label("Median path")
                .legend(move |(x, y)| PathElement::new([(x, y), (x + 20, y)], median_colour.stroke_width(5)));
        }
    }

    let best_fit_colour = RGBColor(255, 20, 147);
    chart
        .draw_series(LineSeries::new(
            grid.best_fit_times
                .iter()
                .copied()
                .zip(grid.best_fit_temperatures.iter().copied()),
            best_fit_colour.stroke_width(4),
        ))
        .map_err(|error| draw_error(output_path, error))?
        .label("Best-fit profile")
        .legend(move |(x, y)| {
            PathElement::new([(x, y), (x + 20, y)], best_fit_colour.stroke_width(4))
        });

    if let Some(profile) = real_profile {
        chart
            .draw_series(DashedLineSeries::new(
                profile
                    .times
                    .iter()
                    .copied()
                    .zip(profile.temperatures.iter().copied()),
                12,
                6,
                BLUE.stroke_width(4),
            ))
            .map_err(|error| draw_error(output_path, error))?
            .label("Known profile")
            .legend(|(x, y)| {
                EmptyElement::at((x, y))
                    + PathElement::new([(0, 0), (8, 0)], BLUE.stroke_width(4))
                    + PathElement::new([(14, 0), (22, 0)], BLUE.stroke_width(4))
            });
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()
        .map_err(|error| draw_error(output_path, error))?;
    draw_probability_scale(&scale_area, output_path, maximum_probability)?;
    root.present()
        .map_err(|error| draw_error(output_path, error))
}

fn probability_colour(intensity: Float) -> RGBColor {
    let intensity = intensity.clamp(0.0, 1.0);
    let stops = [
        (0.0, 0.0, 100.0),     // dark blue
        (0.0, 0.0, 255.0),     // blue
        (173.0, 216.0, 230.0), // light blue
        (0.0, 128.0, 0.0),     // green
        (255.0, 255.0, 0.0),   // yellow
        (255.0, 140.0, 0.0),   // orange
        (255.0, 99.0, 71.0),   // lighter red
        (139.0, 0.0, 0.0),     // dark red
    ];
    let position = intensity * (stops.len() - 1) as Float;
    let index = (position.floor() as usize).min(stops.len() - 2);
    let fraction = position - index as Float;
    let start = stops[index];
    let end = stops[index + 1];
    RGBColor(
        (start.0 + (end.0 - start.0) * fraction).round() as u8,
        (start.1 + (end.1 - start.1) * fraction).round() as u8,
        (start.2 + (end.2 - start.2) * fraction).round() as u8,
    )
}

fn draw_probability_scale<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    output: &Path,
    maximum_probability: Float,
) -> Result<(), PlotError>
where
    DB::ErrorType: std::fmt::Debug,
{
    let mut scale = ChartBuilder::on(area)
        .margin(8)
        .x_label_area_size(0)
        .y_label_area_size(48)
        .build_cartesian_2d(0.0..1.0, 0.0..maximum_probability)
        .map_err(|error| draw_error(output, error))?;
    scale
        .configure_mesh()
        .disable_x_mesh()
        .disable_x_axis()
        .y_labels(5)
        .y_label_formatter(&|value| format!("{:.0}%", value * 100.0))
        .y_desc("Probability")
        .draw()
        .map_err(|error| draw_error(output, error))?;
    scale
        .draw_series((0..100).map(|step| {
            let low = maximum_probability * step as Float / 100.0;
            let high = maximum_probability * (step + 1) as Float / 100.0;
            Rectangle::new(
                [(0.0, low), (1.0, high)],
                probability_colour((step as Float + 0.5) / 100.0).filled(),
            )
        }))
        .map_err(|error| draw_error(output, error))?;
    Ok(())
}

/// Draw all accepted profiles, highlighting frequentist p-values >= 0.5.
///
/// The known profile must use the units stored in the input file.
pub fn plot_accepted_profiles<P>(
    input_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    known_profile: Option<ProfilePath<'_>>,
) -> Result<(), Box<dyn std::error::Error>>
where
    P: serde::de::DeserializeOwned + ProfileCoordinates,
{
    use crate::chrono_outputs::{
        read_all_chrono_batches, AcceptedProfile, AcceptanceInfo,
    };

    let input_path = input_path.as_ref();
    let output_path = output_path.as_ref();

    let reader = read_all_chrono_batches::<AcceptedProfile<P>>(input_path)?;
    let time_unit = reader.time_unit;
    let temp_unit = reader.temp_unit;

    let mut time_min = TimeFloat::INFINITY;
    let mut time_max = TimeFloat::NEG_INFINITY;
    let mut temp_min = Float::INFINITY;
    let mut temp_max = Float::NEG_INFINITY;
    let mut count = 0_u64;

    // First pass: validate profiles and find bounds.
    {
        let mut include = |profile: ProfilePath<'_>| -> Result<(), PlotError> {
            validate_overlay(input_path, profile, time_unit)?;

            for &time in profile.times {
                time_min = time_min.min(time);
                time_max = time_max.max(time);
            }
            for &temperature in profile.temperatures {
                temp_min = temp_min.min(temperature);
                temp_max = temp_max.max(temperature);
            }
            Ok(())
        };

        for batch in reader {
            for record in batch? {
                include(record.profile.coordinates())?;
                count += 1;
            }
        }

        if let Some(profile) = known_profile {
            include(profile)?;
        }
    }

    if count == 0 {
        return Err(invalid(input_path, "no accepted profiles were found").into());
    }

    // Pad the axes, including constant-temperature profiles.
    let time_padding = (time_max - time_min) * 0.02;
    let temp_padding = ((temp_max - temp_min) * 0.05).max(1.0);

    time_min -= time_padding;
    time_max += time_padding;
    temp_min -= temp_padding;
    temp_max += temp_padding;

    let time_range = if time_unit.is_ka_or_ma() {
        time_max..time_min
    } else {
        time_min..time_max
    };

    let root =
        BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("Time-temperature paths ({count} accepted profiles)"),
            ("sans-serif", 28),
        )
        .margin(20)
        .x_label_area_size(55)
        .y_label_area_size(75)
        .build_cartesian_2d(time_range, temp_min..temp_max)?;

    chart
        .configure_mesh()
        .x_desc(format!("Time ({time_unit})"))
        .y_desc(format!("Temperature ({temp_unit})"))
        .draw()?;

    let grey = RGBColor(160, 160, 160);
    let green = RGBColor(0, 150, 0);

    // Draw grey first, then green so highlighted paths remain visible.
    for highlight in [false, true] {
        let reader =
            read_all_chrono_batches::<AcceptedProfile<P>>(input_path)?;
        let mut labelled = false;

        for batch in reader {
            for record in batch? {
                let is_green = matches!(
                    &record.acceptance,
                    AcceptanceInfo::Frequentist { p_test }
                        if *p_test >= 0.5
                );

                if is_green != highlight {
                    continue;
                }

                let profile = record.profile.coordinates();
                let colour = if highlight { green } else { grey };

                // Same segment rendering as plot_time_temperature.
                let series = chart.draw_series(
                    profile.times.windows(2)
                        .zip(profile.temperatures.windows(2))
                        .map(|(times, temperatures)| {
                            PathElement::new(
                                [
                                    (times[0], temperatures[0]),
                                    (times[1], temperatures[1]),
                                ],
                                colour.stroke_width(2),
                            )
                        }),
                )?;

                if !labelled {
                    series
                        .label(if highlight {
                            "Frequentist p >= 0.5"
                        } else {
                            "Other accepted profiles"
                        })
                        .legend(move |(x, y)| {
                            PathElement::new(
                                [(x, y), (x + 20, y)],
                                colour.stroke_width(2),
                            )
                        });
                    labelled = true;
                }
            }
        }
    }

    // Draw the known profile last.
    if let Some(profile) = known_profile {
        chart
            .draw_series(DottedLineSeries::new(
                profile.times.iter().copied()
                    .zip(profile.temperatures.iter().copied()),
                0,
                8,
                |point| Circle::new(point, 2, BLUE.filled()),
            ))?
            .label("Known profile")
            .legend(|(x, y)| {
                EmptyElement::at((x, y))
                    + Circle::new((0, 0), 2, BLUE.filled())
                    + Circle::new((8, 0), 2, BLUE.filled())
                    + Circle::new((16, 0), 2, BLUE.filled())
            });
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.9))
        .border_style(BLACK)
        .draw()?;

    root.present()?;
    Ok(())
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
            "3\n2,0,Ka annum\n0,3,Celsius\n2,0\n0.5,2.5\n0,1\n1,2\n2,0\n",
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
