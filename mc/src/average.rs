// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Streaming consolidation of repeated Monte Carlo trajectories.
//!
//! Fill trajectories are step functions because occupancy changes only at
//! discrete events. Their union of record times is traversed without loading
//! complete files, carrying each repetition's last fill forward before
//! calculating ensemble statistics. Event trajectories are instead grouped
//! into fixed-width bins and normalized into observed frequencies.

use common::charge_transfer::{ElectronicState, Event, RecordedEvent};
use common::numeric::{Float, TimeFloat};
use io::outputs::{
    AverageEventRow, BatchReader, ContinuousValueRow, read_all_batches, write_average_events_csv,
    write_continuous_values_csv,
};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::vec::IntoIter;

/// Default directory containing per-repetition compressed trajectories.
const TEMPORARY_DIRECTORY: &str = "tmp";
/// Prefix used for per-experiment ensemble fill statistics.
const AVERAGE_FILL_PREFIX: &str = "average_fill";
/// Prefix used for per-experiment ensemble event frequencies.
const AVERAGE_EVENT_PREFIX: &str = "average_event";
/// Smallest automatically selected event-count bin width, in seconds.
const MINIMUM_EVENT_BIN_WIDTH: TimeFloat = 0.1;
/// Approximate number of event bins produced when a profile is long enough to
/// require bins wider than [`MINIMUM_EVENT_BIN_WIDTH`].
const TARGET_EVENT_BIN_COUNT: usize = 1_000;
/// Hard limit protecting explicit bin-width requests from excessive output or
/// memory consumption.
const MAXIMUM_EVENT_BIN_COUNT: usize = 1_000_000;

/// Flattens the compressed batches in one result file into individual records.
struct RecordStream {
    /// Source path used when reporting malformed trajectory data.
    path: PathBuf,
    /// Lazy iterator over compressed serialized batches.
    batches: BatchReader<RecordedEvent>,
    /// Remaining records in the currently decoded batch.
    records: IntoIter<RecordedEvent>,
}

impl RecordStream {
    /// Open a trajectory and prepare to decode its first batch lazily.
    fn open(path: PathBuf) -> Result<Self, String> {
        let batches = read_all_batches(&path).map_err(|error| error.to_string())?;
        Ok(Self {
            path,
            batches,
            records: Vec::new().into_iter(),
        })
    }

    /// Return the next record, transparently advancing across batch boundaries.
    fn next_record(&mut self) -> Result<Option<RecordedEvent>, String> {
        loop {
            if let Some(record) = self.records.next() {
                return Ok(Some(record));
            }

            match self.batches.next() {
                Some(Ok(batch)) => self.records = batch.into_iter(),
                Some(Err(error)) => return Err(error.to_string()),
                None => return Ok(None),
            }
        }
    }
}

/// Streaming state for one repetition during a multiway time merge.
struct TrajectoryCursor {
    /// Flattened record source for this repetition.
    stream: RecordStream,
    /// Most recent time incorporated into the ensemble state.
    last_time: TimeFloat,
    /// Fill value carried forward from the most recent record.
    fill: Float,
    /// Look-ahead record used to find the next union time.
    next: Option<RecordedEvent>,
}

/// Monotonic direction shared by all trajectories being averaged.
#[derive(Debug, Clone, Copy)]
enum Direction {
    /// Ordinary experimental time increasing from zero.
    Forward,
    /// Geological age decreasing toward zero.
    Reverse,
}

impl Direction {
    /// Select the earliest next union time in the direction of travel.
    fn next_time(self, current: Option<TimeFloat>, candidate: TimeFloat) -> TimeFloat {
        match current {
            None => candidate,
            Some(current) => match self {
                Self::Forward if candidate.total_cmp(&current).is_lt() => candidate,
                Self::Reverse if candidate.total_cmp(&current).is_gt() => candidate,
                _ => current,
            },
        }
    }

    /// Return whether two consecutive times preserve this direction.
    fn accepts(self, previous: TimeFloat, next: TimeFloat) -> bool {
        match self {
            Self::Forward => next >= previous,
            Self::Reverse => next <= previous,
        }
    }
}

/// Iterator that merges repetitions and emits ensemble fill rows incrementally.
struct AverageFillRows {
    /// One look-ahead cursor per repetition.
    cursors: Vec<TrajectoryCursor>,
    /// Shared monotonic direction inferred from the first changing trajectory.
    direction: Direction,
    /// Initial state, emitted before processing look-ahead records.
    first: Option<ContinuousValueRow>,
    /// Whether all records have been consumed or an error has terminated iteration.
    finished: bool,
}

/// Ensemble summary of fill fractions at one union time.
struct FillStatistics {
    /// Arithmetic mean across repetitions.
    mean: Float,
    /// Population standard deviation across repetitions.
    standard_deviation: Float,
    /// Linearly interpolated 50th percentile.
    median: Float,
    /// Linearly interpolated 10th percentile.
    quantile_0_1: Float,
    /// Linearly interpolated 90th percentile.
    quantile_0_9: Float,
    /// Linearly interpolated lower quartile.
    quantile_0_25: Float,
    /// Linearly interpolated upper quartile.
    quantile_0_75: Float,
}

impl FillStatistics {
    /// Sort fill values and calculate population statistics.
    fn from_values(mut values: Vec<Float>) -> Self {
        debug_assert!(!values.is_empty());
        values.sort_by(Float::total_cmp);

        let mean = values.iter().sum::<Float>() / values.len() as Float;
        let variance = values
            .iter()
            .map(|value| (value - mean).powi(2))
            .sum::<Float>()
            / values.len() as Float;

        Self {
            mean,
            standard_deviation: variance.sqrt(),
            median: quantile(&values, 0.5),
            quantile_0_1: quantile(&values, 0.1),
            quantile_0_9: quantile(&values, 0.9),
            quantile_0_25: quantile(&values, 0.25),
            quantile_0_75: quantile(&values, 0.75),
        }
    }

    /// Combine these fill statistics with their physical time and temperature.
    fn row(self, time: TimeFloat, temperature: Float) -> ContinuousValueRow {
        ContinuousValueRow {
            time,
            temperature,
            fill: self.mean,
            fill_standard_deviation: self.standard_deviation,
            fill_median: self.median,
            fill_quantile_0_1: self.quantile_0_1,
            fill_quantile_0_9: self.quantile_0_9,
            fill_quantile_0_25: self.quantile_0_25,
            fill_quantile_0_75: self.quantile_0_75,
        }
    }
}

/// Calculate a linearly interpolated quantile from a non-empty sorted slice.
fn quantile(sorted_values: &[Float], probability: Float) -> Float {
    let position = probability * (sorted_values.len() - 1) as Float;
    let lower_index = position.floor() as usize;
    let upper_index = position.ceil() as usize;
    let fraction = position - lower_index as Float;
    sorted_values[lower_index]
        + fraction * (sorted_values[upper_index] - sorted_values[lower_index])
}

impl AverageFillRows {
    /// Open all trajectories, validate their initial states, and infer direction.
    ///
    /// Every trajectory must begin at exactly the same time. Temperature is
    /// averaged only over records present at a union time, whereas fill is
    /// carried forward independently for every repetition.
    fn new(paths: Vec<PathBuf>) -> Result<Self, String> {
        let mut cursors = Vec::with_capacity(paths.len());
        let mut initial_time = None;
        let mut initial_temperature_sum = 0.0;
        let mut initial_fills = Vec::with_capacity(paths.len());

        for path in paths {
            let mut stream = RecordStream::open(path.clone())?;
            let first = stream
                .next_record()?
                .ok_or_else(|| format!("temporary result file {} is empty", path.display()))?;
            validate_record(&path, &first)?;

            if let Some(expected) = initial_time {
                if first.time != expected {
                    return Err(format!(
                        "temporary result files do not share one initial time: {} starts at {}, expected {expected}",
                        path.display(),
                        first.time,
                    ));
                }
            } else {
                initial_time = Some(first.time);
            }

            initial_temperature_sum += first.temperature;
            initial_fills.push(first.fill);
            let next = stream.next_record()?;
            if let Some(record) = &next {
                validate_record(&path, record)?;
            }
            cursors.push(TrajectoryCursor {
                stream,
                last_time: first.time,
                fill: first.fill,
                next,
            });
        }

        let direction = cursors
            .iter()
            .filter_map(|cursor| {
                let next_time = cursor.next.as_ref()?.time;
                match next_time.total_cmp(&cursor.last_time) {
                    Ordering::Greater => Some(Direction::Forward),
                    Ordering::Less => Some(Direction::Reverse),
                    Ordering::Equal => None,
                }
            })
            .next()
            .unwrap_or(Direction::Forward);

        for cursor in &cursors {
            if let Some(next) = &cursor.next
                && !direction.accepts(cursor.last_time, next.time)
            {
                return Err(format!(
                    "times in {} do not follow the same direction as the other repetitions",
                    cursor.stream.path.display()
                ));
            }
        }

        let count = cursors.len() as Float;
        let first = FillStatistics::from_values(initial_fills).row(
            initial_time.expect("at least one path is required"),
            initial_temperature_sum / count,
        );

        Ok(Self {
            cursors,
            direction,
            first: Some(first),
            finished: false,
        })
    }

    /// Consume all consecutive records for one cursor at an exact union time.
    ///
    /// The last fill at that time becomes the carried state; temperatures from
    /// all same-time records contribute to the returned sum and count.
    fn advance_cursor_at(
        cursor: &mut TrajectoryCursor,
        time: TimeFloat,
        direction: Direction,
    ) -> Result<(Float, usize), String> {
        let mut temperature_sum = 0.0;
        let mut temperature_count = 0;

        while cursor
            .next
            .as_ref()
            .is_some_and(|record| record.time == time)
        {
            let record = cursor.next.take().expect("record presence was checked");
            cursor.last_time = record.time;
            cursor.fill = record.fill;
            temperature_sum += record.temperature;
            temperature_count += 1;
            cursor.next = cursor.stream.next_record()?;

            if let Some(next) = &cursor.next {
                validate_record(&cursor.stream.path, next)?;
                if !direction.accepts(cursor.last_time, next.time) {
                    return Err(format!(
                        "times in {} are not monotonic",
                        cursor.stream.path.display()
                    ));
                }
            }
        }

        Ok((temperature_sum, temperature_count))
    }
}

impl Iterator for AverageFillRows {
    type Item = Result<ContinuousValueRow, String>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        if let Some(first) = self.first.take() {
            return Some(Ok(first));
        }

        let mut time = None;
        for cursor in &self.cursors {
            if let Some(record) = &cursor.next {
                time = Some(self.direction.next_time(time, record.time));
            }
        }
        let Some(time) = time else {
            self.finished = true;
            return None;
        };

        let mut temperature_sum = 0.0;
        let mut temperature_count = 0usize;
        for cursor in &mut self.cursors {
            match Self::advance_cursor_at(cursor, time, self.direction) {
                Ok((cursor_temperature_sum, cursor_temperature_count)) => {
                    temperature_sum += cursor_temperature_sum;
                    temperature_count += cursor_temperature_count;
                }
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            }
        }
        if temperature_count == 0 {
            self.finished = true;
            return Some(Err(format!("no temperature was recorded at time {time}")));
        }

        let fills = self.cursors.iter().map(|cursor| cursor.fill).collect();
        Some(Ok(FillStatistics::from_values(fills)
            .row(time, temperature_sum / temperature_count as Float)))
    }
}

/// Raw event counters for one fixed interval of simulation time.
///
/// Events at a positive boundary are assigned to the bin ending at that
/// boundary. Counters are converted to per-second, per-repetition frequencies
/// only after every trajectory has been read.
struct EventBin {
    /// Inclusive left edge of the bin, in seconds.
    start_time: TimeFloat,
    /// Right edge and output timestamp of the bin, in seconds.
    end_time: TimeFloat,

    /// Ground-state localised recombination events.
    localised_recombination_ground_count: usize,
    /// Excited-state localised recombination events.
    localised_recombination_excited_count: usize,
    /// Ground-state conduction-band recombination events.
    delocalised_recombination_ground_count: usize,
    /// Excited-state conduction-band recombination events.
    delocalised_recombination_excited_count: usize,

    /// Ground-state localised trap-to-trap events.
    localised_retrapping_ground_count: usize,
    /// Excited-state localised trap-to-trap events.
    localised_retrapping_excited_count: usize,
    /// Ground-state conduction-band retrapping events.
    delocalised_retrapping_ground_count: usize,
    /// Excited-state conduction-band retrapping events.
    delocalised_retrapping_excited_count: usize,

    /// Total events whose electron originated from a ground state.
    ground_count: usize,
    /// Total events whose electron originated from an excited state.
    excited_count: usize,

    /// Total localised and delocalised recombination events.
    recombination_count: usize,
    /// Total localised and delocalised retrapping events.
    retrapping_count: usize,

    /// Total irradiation-driven filling events.
    filling_count: usize,
}

impl EventBin {
    /// Create an empty bin at `origin + index * width`.
    fn new(index: usize, origin: TimeFloat, width: TimeFloat) -> Self {
        Self {
            start_time: origin + index as TimeFloat * width,
            end_time: origin + (index + 1) as TimeFloat * width,
            localised_recombination_ground_count: 0,
            localised_recombination_excited_count: 0,
            delocalised_recombination_ground_count: 0,
            delocalised_recombination_excited_count: 0,
            localised_retrapping_ground_count: 0,
            localised_retrapping_excited_count: 0,
            delocalised_retrapping_ground_count: 0,
            delocalised_retrapping_excited_count: 0,
            ground_count: 0,
            excited_count: 0,
            recombination_count: 0,
            retrapping_count: 0,
            filling_count: 0,
        }
    }

    /// Increment the specific pathway and its aggregate state/category counters.
    fn record(&mut self, event: Event) {
        match event {
            Event::LocalisedRecombination {
                state: ElectronicState::Ground,
                ..
            } => {
                self.localised_recombination_ground_count += 1;
                self.ground_count += 1;
                self.recombination_count += 1;
            }
            Event::LocalisedRecombination {
                state: ElectronicState::Excited,
                ..
            } => {
                self.localised_recombination_excited_count += 1;
                self.excited_count += 1;
                self.recombination_count += 1;
            }
            Event::DelocalisedRecombination {
                state: ElectronicState::Ground,
                ..
            } => {
                self.delocalised_recombination_ground_count += 1;
                self.ground_count += 1;
                self.recombination_count += 1;
            }
            Event::DelocalisedRecombination {
                state: ElectronicState::Excited,
                ..
            } => {
                self.delocalised_recombination_excited_count += 1;
                self.excited_count += 1;
                self.recombination_count += 1;
            }
            Event::LocalisedRetrapping {
                state: ElectronicState::Ground,
                ..
            } => {
                self.localised_retrapping_ground_count += 1;
                self.ground_count += 1;
                self.retrapping_count += 1;
            }
            Event::LocalisedRetrapping {
                state: ElectronicState::Excited,
                ..
            } => {
                self.localised_retrapping_excited_count += 1;
                self.excited_count += 1;
                self.retrapping_count += 1;
            }
            Event::DelocalisedRetrapping {
                state: ElectronicState::Ground,
                ..
            } => {
                self.delocalised_retrapping_ground_count += 1;
                self.ground_count += 1;
                self.retrapping_count += 1;
            }
            Event::DelocalisedRetrapping {
                state: ElectronicState::Excited,
                ..
            } => {
                self.delocalised_retrapping_excited_count += 1;
                self.excited_count += 1;
                self.retrapping_count += 1;
            }

            Event::FillingStandard { .. }
            | Event::FillingLoss { .. }
            | Event::FillingHoleOnly { .. }
            | Event::FillingTrapOnly { .. } => self.filling_count += 1,

            Event::None => {}
            _ => {}
        }
    }

    /// Normalize raw counts by repetition count.
    fn averaged(self, repetition_count: usize) -> AverageEventRow {
        debug_assert!(self.end_time > self.start_time);
        let repetitions = repetition_count as Float;
        AverageEventRow {
            time: self.end_time,
            localised_recombination_ground_count: self.localised_recombination_ground_count
                as Float
                / repetitions,
            localised_recombination_excited_count: self.localised_recombination_excited_count
                as Float
                / repetitions,
            delocalised_recombination_ground_count: self.delocalised_recombination_ground_count
                as Float
                / repetitions,
            delocalised_recombination_excited_count: self.delocalised_recombination_excited_count
                as Float
                / repetitions,
            localised_retrapping_ground_count: self.localised_retrapping_ground_count as Float
                / repetitions,
            localised_retrapping_excited_count: self.localised_retrapping_excited_count as Float
                / repetitions,
            delocalised_retrapping_ground_count: self.delocalised_retrapping_ground_count as Float
                / repetitions,
            delocalised_retrapping_excited_count: self.delocalised_retrapping_excited_count
                as Float
                / repetitions,
            ground_count: self.ground_count as Float / repetitions,
            excited_count: self.excited_count as Float / repetitions,
            recombination_count: self.recombination_count as Float / repetitions,
            retrapping_count: self.retrapping_count as Float / repetitions,
            filling_count: self.filling_count as Float / repetitions,
        }
    }
}

/// Validate the coordinates needed for event binning.
fn validate_event_record(path: &Path, record: &RecordedEvent) -> Result<(), String> {
    validate_record(path, record)?;
    if record.time < 0.0 {
        return Err(format!(
            "{} contains a negative time {}, which cannot be placed in non-negative bins",
            path.display(),
            record.time,
        ));
    }
    Ok(())
}

/// Find the complete time extent represented by a set of trajectories.
fn event_time_bounds(paths: &[PathBuf]) -> Result<(TimeFloat, TimeFloat), String> {
    let mut minimum = TimeFloat::INFINITY;
    let mut maximum = TimeFloat::NEG_INFINITY;

    for path in paths {
        let batches = read_all_batches::<RecordedEvent>(path).map_err(|error| error.to_string())?;
        for batch in batches {
            for record in batch.map_err(|error| error.to_string())? {
                validate_event_record(path, &record)?;
                minimum = minimum.min(record.time);
                maximum = maximum.max(record.time);
            }
        }
    }

    if !minimum.is_finite() {
        return Err("temporary event trajectories contain no records".to_string());
    }
    Ok((minimum, maximum))
}

/// Select a requested width or an automatic width bounded to roughly one
/// thousand bins. Short simulations retain the historical 0.1-second bins.
fn event_bin_width(
    minimum_time: TimeFloat,
    maximum_time: TimeFloat,
    requested_width: Option<TimeFloat>,
) -> Result<TimeFloat, String> {
    if let Some(width) = requested_width {
        if !width.is_finite() || width <= 0.0 {
            return Err(format!(
                "event-bin width must be finite and positive, got {width}"
            ));
        }
        return Ok(width);
    }

    let duration = maximum_time - minimum_time;
    Ok(MINIMUM_EVENT_BIN_WIDTH.max(duration / TARGET_EVENT_BIN_COUNT as TimeFloat))
}

/// Calculate and validate the number of bins before attempting an allocation.
fn event_bin_count(
    minimum_time: TimeFloat,
    maximum_time: TimeFloat,
    width: TimeFloat,
) -> Result<usize, String> {
    let count = ((maximum_time - minimum_time) / width).ceil().max(1.0);
    if !count.is_finite() || count > MAXIMUM_EVENT_BIN_COUNT as TimeFloat {
        return Err(format!(
            "event-bin width {width} seconds would create {count:.0} bins; the maximum is {MAXIMUM_EVENT_BIN_COUNT}"
        ));
    }
    Ok(count as usize)
}

/// Return the all-zero event-count row at the left edge of the first bin.
fn initial_event_row(time: TimeFloat) -> AverageEventRow {
    AverageEventRow {
        time,
        localised_recombination_ground_count: 0.0,
        localised_recombination_excited_count: 0.0,
        delocalised_recombination_ground_count: 0.0,
        delocalised_recombination_excited_count: 0.0,
        localised_retrapping_ground_count: 0.0,
        localised_retrapping_excited_count: 0.0,
        delocalised_retrapping_ground_count: 0.0,
        delocalised_retrapping_excited_count: 0.0,
        ground_count: 0.0,
        excited_count: 0.0,
        recombination_count: 0.0,
        retrapping_count: 0.0,
        filling_count: 0.0,
    }
}
/// Validate finite physical coordinates and a bounded fill fraction.
fn validate_record(path: &Path, record: &RecordedEvent) -> Result<(), String> {
    if !record.time.is_finite() {
        return Err(format!("{} contains a non-finite time", path.display()));
    }
    if !record.temperature.is_finite() {
        return Err(format!(
            "{} contains a non-finite temperature",
            path.display()
        ));
    }
    if !record.fill.is_finite() || !(0.0..=1.0).contains(&record.fill) {
        return Err(format!(
            "{} contains an invalid fill value {}",
            path.display(),
            record.fill
        ));
    }
    Ok(())
}

/// Discover result files and group them by experiment index.
fn temporary_result_paths(directory: &Path) -> Result<BTreeMap<usize, Vec<PathBuf>>, String> {
    let entries = fs::read_dir(directory).map_err(|error| {
        format!(
            "failed to read temporary result directory {}: {error}",
            directory.display()
        )
    })?;
    let mut paths_by_experiment = BTreeMap::<usize, Vec<(usize, PathBuf)>>::new();

    for entry in entries {
        let entry = entry.map_err(|error| {
            format!(
                "failed to read an entry in {}: {error}",
                directory.display()
            )
        })?;
        let path = entry.path();
        let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };
        if !path.is_file() {
            continue;
        }
        let Some(indices) = name
            .strip_prefix("experiment_results_")
            .and_then(|name| name.strip_suffix(".bin.gz"))
        else {
            continue;
        };
        let (experiment_index, repetition_index) = indices.split_once('_').ok_or_else(|| {
            format!(
                "temporary result file {} does not contain experiment and repetition indices",
                path.display()
            )
        })?;
        let experiment_index = experiment_index.parse::<usize>().map_err(|error| {
            format!(
                "temporary result file {} has an invalid experiment index: {error}",
                path.display()
            )
        })?;
        let repetition_index = repetition_index.parse::<usize>().map_err(|error| {
            format!(
                "temporary result file {} has an invalid repetition index: {error}",
                path.display()
            )
        })?;
        paths_by_experiment
            .entry(experiment_index)
            .or_default()
            .push((repetition_index, path));
    }

    if paths_by_experiment.is_empty() {
        return Err(format!(
            "no Monte Carlo temporary result files were found in {}",
            directory.display()
        ));
    }

    Ok(paths_by_experiment
        .into_iter()
        .map(|(experiment_index, mut paths)| {
            paths.sort_by_key(|(repetition_index, _)| *repetition_index);
            (
                experiment_index,
                paths.into_iter().map(|(_, path)| path).collect(),
            )
        })
        .collect())
}

/// Construct a per-experiment CSV path in the requested output directory.
fn average_output_path(directory: &Path, prefix: &str, experiment_index: usize) -> PathBuf {
    directory.join(format!("{prefix}_{experiment_index}.csv"))
}

/// Average one experiment's continuous fill trajectories into one CSV file.
fn average_fill_paths(paths: Vec<PathBuf>, output_file: &Path) -> Result<(), String> {
    let rows = AverageFillRows::new(paths)?;
    write_continuous_values_csv(output_file, rows).map_err(|error| error.to_string())
}

/// Average one experiment's event trajectories into one CSV file.
fn average_event_paths(
    paths: Vec<PathBuf>,
    output_file: &Path,
    requested_bin_width: Option<TimeFloat>,
) -> Result<(), String> {
    let repetition_count = paths.len();
    let (minimum_time, maximum_time) = event_time_bounds(&paths)?;
    let bin_width = event_bin_width(minimum_time, maximum_time, requested_bin_width)?;
    let bin_count = event_bin_count(minimum_time, maximum_time, bin_width)?;
    let mut bins = Vec::<EventBin>::new();
    bins.try_reserve_exact(bin_count)
        .map_err(|error| format!("failed to allocate {bin_count} event bins: {error}"))?;
    for index in 0..bin_count {
        bins.push(EventBin::new(index, minimum_time, bin_width));
    }

    for path in paths {
        let batches =
            read_all_batches::<RecordedEvent>(&path).map_err(|error| error.to_string())?;
        for batch in batches {
            for record in batch.map_err(|error| error.to_string())? {
                validate_event_record(&path, &record)?;
                if record.event == Event::None {
                    continue;
                }

                let relative_time = record.time - minimum_time;
                let endpoint = (relative_time / bin_width).ceil() as usize;
                let bin_index = endpoint.saturating_sub(1).min(bin_count - 1);
                bins[bin_index].record(record.event);
            }
        }
    }

    let rows = std::iter::once(initial_event_row(minimum_time))
        .chain(bins.into_iter().map(|bin| bin.averaged(repetition_count)))
        .map(Ok::<_, std::convert::Infallible>);
    write_average_events_csv(output_file, rows).map_err(|error| error.to_string())
}

/// Average the continuous fill trajectories in `tmp/` and write one
/// `average_fill_<experiment>.csv` file per experiment.
pub fn average_fill() -> Result<(), String> {
    average_fill_in(TEMPORARY_DIRECTORY, ".")
}

/// Average temporary trajectories from `temporary_directory` into one CSV per
/// experiment in `output_directory`.
///
/// Output times are the union of all repetition timestamps. A repetition with
/// no event at a particular union time contributes its most recently observed
/// fill, reflecting piecewise-constant occupancy between events.
pub fn average_fill_in(
    temporary_directory: impl AsRef<Path>,
    output_directory: impl AsRef<Path>,
) -> Result<(), String> {
    for (experiment_index, paths) in temporary_result_paths(temporary_directory.as_ref())? {
        let output_file = average_output_path(
            output_directory.as_ref(),
            AVERAGE_FILL_PREFIX,
            experiment_index,
        );
        average_fill_paths(paths, &output_file)?;
    }
    Ok(())
}

/// Average event counts from `tmp/` into duration-aware bins and write one
/// `average_event_<experiment>.csv` file per experiment.
///
/// Short runs retain 0.1-second bins. Longer runs target approximately one
/// thousand bins.
pub fn average_events() -> Result<(), String> {
    average_events_in(TEMPORARY_DIRECTORY, ".")
}

/// Average event counts using an explicit bin width in seconds.
pub fn average_events_with_bin_width(bin_width: TimeFloat) -> Result<(), String> {
    average_events_in_with_bin_width(TEMPORARY_DIRECTORY, ".", bin_width)
}

/// Average temporary event trajectories from `temporary_directory` into one
/// CSV per experiment in `output_directory`.
///
/// Counts are divided by the number of input trajectories. The plotting layer
/// divides these averaged bin counts by the applicable bin width to produce
/// observed frequencies. Short runs retain 0.1-second bins, while longer runs
/// use a bounded duration-aware width.
pub fn average_events_in(
    temporary_directory: impl AsRef<Path>,
    output_directory: impl AsRef<Path>,
) -> Result<(), String> {
    average_events_in_configured(temporary_directory, output_directory, None)
}

/// Average temporary event trajectories using an explicit bin width in
/// seconds.
pub fn average_events_in_with_bin_width(
    temporary_directory: impl AsRef<Path>,
    output_directory: impl AsRef<Path>,
    bin_width: TimeFloat,
) -> Result<(), String> {
    average_events_in_configured(temporary_directory, output_directory, Some(bin_width))
}

/// Shared implementation for automatic and explicitly sized event bins.
fn average_events_in_configured(
    temporary_directory: impl AsRef<Path>,
    output_directory: impl AsRef<Path>,
    requested_bin_width: Option<TimeFloat>,
) -> Result<(), String> {
    for (experiment_index, paths) in temporary_result_paths(temporary_directory.as_ref())? {
        let output_file = average_output_path(
            output_directory.as_ref(),
            AVERAGE_EVENT_PREFIX,
            experiment_index,
        );
        average_event_paths(paths, &output_file, requested_bin_width)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::charge_transfer::Event;
    use common::constants::time::{self, TimeUnit};
    use common::place_ids::PlaceId;
    use io::outputs::{
        append_monte_carlo_experiment_batch_to_file, create_monte_carlo_experiment_file,
    };
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
    use std::time::{SystemTime, UNIX_EPOCH};

    static NEXT_TEMPORARY_DIRECTORY: AtomicU64 = AtomicU64::new(0);

    fn temporary_directory() -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        let sequence = NEXT_TEMPORARY_DIRECTORY.fetch_add(1, AtomicOrdering::Relaxed);
        let directory = std::env::temp_dir().join(format!(
            "mcrustlum_average_{}_{}_{}",
            std::process::id(),
            unique,
            sequence,
        ));
        fs::create_dir(&directory).unwrap();
        directory
    }

    fn record(time: TimeFloat, temperature: Float, fill: Float) -> RecordedEvent {
        RecordedEvent {
            time,
            temperature,
            fill,
            event: Event::None,
        }
    }

    fn event_record(time: TimeFloat, event: Event) -> RecordedEvent {
        RecordedEvent {
            time,
            temperature: 300.0,
            fill: 0.5,
            event,
        }
    }

    fn write_records(path: &Path, records: &[RecordedEvent]) {
        create_monte_carlo_experiment_file(path).unwrap();
        for batch in records.chunks(2) {
            append_monte_carlo_experiment_batch_to_file(path, batch).unwrap();
        }
    }

    #[test]
    fn unions_times_and_carries_each_fill_forward_before_averaging() {
        let directory = temporary_directory();
        let temporary_directory = directory.join("tmp");
        let output_file = directory.join("average_fill_0.csv");
        fs::create_dir(&temporary_directory).unwrap();

        write_records(
            &temporary_directory.join("experiment_results_0_0.bin.gz"),
            &[
                record(0.0, 100.0, 0.2),
                record(2.0, 120.0, 0.4),
                record(4.0, 140.0, 0.6),
            ],
        );
        write_records(
            &temporary_directory.join("experiment_results_0_1.bin.gz"),
            &[
                record(0.0, 100.0, 0.6),
                record(1.0, 110.0, 0.8),
                record(3.0, 130.0, 0.2),
                record(4.0, 140.0, 0.4),
            ],
        );

        average_fill_in(&temporary_directory, &directory).unwrap();
        let contents = fs::read_to_string(&output_file).unwrap();
        fs::remove_dir_all(directory).unwrap();

        let mut lines = contents.lines();
        assert_eq!(
            lines.next().unwrap(),
            "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9,fill_quantile_0_25,fill_quantile_0_75"
        );
        let rows = lines
            .map(|line| {
                line.split(',')
                    .map(|value| value.parse::<Float>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let expected = [
            [0.0, 100.0, 0.4, 0.2, 0.4, 0.24, 0.56, 0.3, 0.5],
            [1.0, 110.0, 0.5, 0.3, 0.5, 0.26, 0.74, 0.35, 0.65],
            [2.0, 120.0, 0.6, 0.2, 0.6, 0.44, 0.76, 0.5, 0.7],
            [3.0, 130.0, 0.3, 0.1, 0.3, 0.22, 0.38, 0.25, 0.35],
            [4.0, 140.0, 0.5, 0.1, 0.5, 0.42, 0.58, 0.45, 0.55],
        ];

        assert_eq!(rows.len(), expected.len());
        for (row, expected_row) in rows.iter().zip(expected) {
            assert_eq!(row.len(), expected_row.len());
            for (actual, expected) in row.iter().zip(expected_row) {
                assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
            }
        }
    }

    #[test]
    fn preserves_reverse_simulation_order() {
        let directory = temporary_directory();
        let temporary_directory = directory.join("tmp");
        let output_file = directory.join("average_fill_0.csv");
        fs::create_dir(&temporary_directory).unwrap();

        write_records(
            &temporary_directory.join("experiment_results_0_0.bin.gz"),
            &[
                record(4.0, 140.0, 0.2),
                record(2.0, 120.0, 0.4),
                record(0.0, 100.0, 0.6),
            ],
        );
        write_records(
            &temporary_directory.join("experiment_results_0_1.bin.gz"),
            &[
                record(4.0, 140.0, 0.6),
                record(3.0, 130.0, 0.8),
                record(1.0, 110.0, 0.2),
                record(0.0, 100.0, 0.4),
            ],
        );

        average_fill_in(&temporary_directory, &directory).unwrap();
        let times = fs::read_to_string(&output_file)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                line.split(',')
                    .next()
                    .unwrap()
                    .parse::<TimeFloat>()
                    .unwrap()
            })
            .collect::<Vec<_>>();
        fs::remove_dir_all(directory).unwrap();

        assert_eq!(times, vec![4.0, 3.0, 2.0, 1.0, 0.0]);
    }

    #[test]
    fn bins_events_and_averages_counts_over_repetitions() {
        let directory = temporary_directory();
        let temporary_directory = directory.join("tmp");
        let output_file = directory.join("average_event_0.csv");
        fs::create_dir(&temporary_directory).unwrap();
        let place = PlaceId::new(0).unwrap();

        write_records(
            &temporary_directory.join("experiment_results_0_0.bin.gz"),
            &[
                event_record(0.0, Event::None),
                event_record(
                    0.0005,
                    Event::LocalisedRecombination {
                        source: place,
                        hole: place,
                        state: ElectronicState::Ground,
                    },
                ),
                event_record(
                    0.001,
                    Event::DelocalisedRecombination {
                        source: place,
                        hole: place,
                        state: ElectronicState::Ground,
                    },
                ),
                event_record(
                    0.0015,
                    Event::FillingStandard {
                        trap: place,
                        hole: place,
                    },
                ),
                event_record(0.003, Event::None),
            ],
        );
        write_records(
            &temporary_directory.join("experiment_results_0_1.bin.gz"),
            &[
                event_record(0.0, Event::None),
                event_record(
                    0.0002,
                    Event::LocalisedRecombination {
                        source: place,
                        hole: place,
                        state: ElectronicState::Ground,
                    },
                ),
                event_record(
                    0.0004,
                    Event::LocalisedRecombination {
                        source: place,
                        hole: place,
                        state: ElectronicState::Ground,
                    },
                ),
                event_record(
                    0.0012,
                    Event::DelocalisedRecombination {
                        source: place,
                        hole: place,
                        state: ElectronicState::Excited,
                    },
                ),
                event_record(0.003, Event::None),
            ],
        );

        average_events_in(&temporary_directory, &directory).unwrap();
        let contents = fs::read_to_string(&output_file).unwrap();
        fs::remove_dir_all(directory).unwrap();
        let rows = contents
            .lines()
            .skip(1)
            .map(|line| line.split(',').collect::<Vec<_>>())
            .collect::<Vec<_>>();

        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0], vec!["0"; 14]);
        assert_eq!(rows[1][0], "0.1");
        assert_eq!(rows[1][1], "1.5");
        assert_eq!(rows[1][3], "0.5");
        assert_eq!(rows[1][4], "0.5");
        assert_eq!(rows[1][9], "2");
        assert_eq!(rows[1][10], "0.5");
        assert_eq!(rows[1][11], "2.5");
        assert_eq!(rows[1][13], "0.5");
    }

    #[test]
    fn bounds_automatic_bins_for_geological_times() {
        let directory = temporary_directory();
        let temporary_directory = directory.join("tmp");
        let output_file = directory.join("average_event_0.csv");
        fs::create_dir(&temporary_directory).unwrap();
        let place = PlaceId::new(0).unwrap();
        let half_ma = time::convert_to_seconds(TimeUnit::MaAnnum, 0.5).unwrap();

        write_records(
            &temporary_directory.join("experiment_results_0_0.bin.gz"),
            &[
                event_record(half_ma, Event::None),
                event_record(
                    half_ma / 2.0,
                    Event::FillingStandard {
                        trap: place,
                        hole: place,
                    },
                ),
                event_record(0.0, Event::None),
            ],
        );

        average_events_in(&temporary_directory, &directory).unwrap();
        let contents = fs::read_to_string(&output_file).unwrap();
        fs::remove_dir_all(directory).unwrap();
        let rows = contents.lines().skip(1).collect::<Vec<_>>();

        assert_eq!(rows.len(), TARGET_EVENT_BIN_COUNT + 1);
        let final_time = rows
            .last()
            .unwrap()
            .split(',')
            .next()
            .unwrap()
            .parse::<TimeFloat>()
            .unwrap();
        assert!((final_time - half_ma).abs() < 1.0);
        let filling_total = rows
            .iter()
            .map(|row| row.split(',').nth(13).unwrap().parse::<Float>().unwrap())
            .sum::<Float>();
        assert_eq!(filling_total, 1.0);
    }

    #[test]
    fn rejects_an_explicit_width_that_would_create_too_many_bins() {
        let half_ma = time::convert_to_seconds(TimeUnit::MaAnnum, 0.5).unwrap();
        let error = event_bin_count(0.0, half_ma, MINIMUM_EVENT_BIN_WIDTH).unwrap_err();

        assert!(error.contains("the maximum is 1000000"));
    }

    #[test]
    fn separates_fill_and_event_outputs_by_experiment() {
        let directory = temporary_directory();
        let temporary_directory = directory.join("tmp");
        fs::create_dir(&temporary_directory).unwrap();
        let place = PlaceId::new(0).unwrap();

        for repetition_index in 0..2 {
            write_records(
                &temporary_directory
                    .join(format!("experiment_results_0_{repetition_index}.bin.gz")),
                &[
                    event_record(0.0, Event::None),
                    RecordedEvent {
                        time: 0.05,
                        temperature: 300.0,
                        fill: 0.2,
                        event: Event::FillingStandard {
                            trap: place,
                            hole: place,
                        },
                    },
                ],
            );
            write_records(
                &temporary_directory
                    .join(format!("experiment_results_1_{repetition_index}.bin.gz")),
                &[
                    event_record(0.0, Event::None),
                    RecordedEvent {
                        time: 0.05,
                        temperature: 300.0,
                        fill: 0.8,
                        event: Event::None,
                    },
                ],
            );
        }

        average_fill_in(&temporary_directory, &directory).unwrap();
        average_events_in(&temporary_directory, &directory).unwrap();

        let last_fill = |experiment_index| {
            fs::read_to_string(directory.join(format!("average_fill_{experiment_index}.csv")))
                .unwrap()
                .lines()
                .last()
                .unwrap()
                .split(',')
                .nth(2)
                .unwrap()
                .parse::<Float>()
                .unwrap()
        };
        let filling_count = |experiment_index| {
            fs::read_to_string(directory.join(format!("average_event_{experiment_index}.csv")))
                .unwrap()
                .lines()
                .last()
                .unwrap()
                .split(',')
                .nth(13)
                .unwrap()
                .parse::<Float>()
                .unwrap()
        };

        assert_eq!(last_fill(0), 0.2);
        assert_eq!(last_fill(1), 0.8);
        assert_eq!(filling_count(0), 1.0);
        assert_eq!(filling_count(1), 0.0);
        fs::remove_dir_all(directory).unwrap();
    }
}
