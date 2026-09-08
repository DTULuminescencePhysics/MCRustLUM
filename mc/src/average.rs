// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Consolidation of repeated Monte Carlo trajectories.

use common::charge_transfer::{ElectronicState, Event, RecordedEvent};
use common::numeric::{Float, TimeFloat};
use io::outputs::{
    AverageEventRow, BatchReader, ContinuousValueRow, read_all_batches, write_average_events_csv,
    write_continuous_values_csv,
};
use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};
use std::vec::IntoIter;

const TEMPORARY_DIRECTORY: &str = "tmp";
const AVERAGE_FILL_FILE: &str = "average_fill.csv";
const AVERAGE_EVENT_FILE: &str = "average_event.csv";
const EVENT_BIN_WIDTH: TimeFloat = 0.1;

struct RecordStream {
    path: PathBuf,
    batches: BatchReader<RecordedEvent>,
    records: IntoIter<RecordedEvent>,
}

impl RecordStream {
    fn open(path: PathBuf) -> Result<Self, String> {
        let batches = read_all_batches(&path).map_err(|error| error.to_string())?;
        Ok(Self {
            path,
            batches,
            records: Vec::new().into_iter(),
        })
    }

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

struct TrajectoryCursor {
    stream: RecordStream,
    last_time: TimeFloat,
    fill: Float,
    next: Option<RecordedEvent>,
}

#[derive(Debug, Clone, Copy)]
enum Direction {
    Forward,
    Reverse,
}

impl Direction {
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

    fn accepts(self, previous: TimeFloat, next: TimeFloat) -> bool {
        match self {
            Self::Forward => next >= previous,
            Self::Reverse => next <= previous,
        }
    }
}

struct AverageFillRows {
    cursors: Vec<TrajectoryCursor>,
    direction: Direction,
    first: Option<ContinuousValueRow>,
    finished: bool,
}

struct FillStatistics {
    mean: Float,
    standard_deviation: Float,
    median: Float,
    quantile_0_1: Float,
    quantile_0_9: Float,
}

impl FillStatistics {
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
        }
    }

    fn row(self, time: TimeFloat, temperature: Float) -> ContinuousValueRow {
        ContinuousValueRow {
            time,
            temperature,
            fill: self.mean,
            fill_standard_deviation: self.standard_deviation,
            fill_median: self.median,
            fill_quantile_0_1: self.quantile_0_1,
            fill_quantile_0_9: self.quantile_0_9,
        }
    }
}

fn quantile(sorted_values: &[Float], probability: Float) -> Float {
    let position = probability * (sorted_values.len() - 1) as Float;
    let lower_index = position.floor() as usize;
    let upper_index = position.ceil() as usize;
    let fraction = position - lower_index as Float;
    sorted_values[lower_index]
        + fraction * (sorted_values[upper_index] - sorted_values[lower_index])
}

impl AverageFillRows {
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

struct EventBin {
    start_time: TimeFloat,
    end_time: TimeFloat,

    localised_recombination_ground_count: usize,
    localised_recombination_excited_count: usize,
    delocalised_recombination_ground_count: usize,
    delocalised_recombination_excited_count: usize,

    localised_retrapping_ground_count: usize,
    localised_retrapping_excited_count: usize,
    delocalised_retrapping_ground_count: usize,
    delocalised_retrapping_excited_count: usize,

    ground_count: usize,
    excited_count: usize,

    recombination_count: usize,
    retrapping_count: usize,

    filling_count: usize,
}

impl EventBin {
    fn new(index: usize) -> Self {
        Self {
            start_time: index as TimeFloat * EVENT_BIN_WIDTH,
            end_time: (index + 1) as TimeFloat * EVENT_BIN_WIDTH,
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

            Event::Filling { .. } => self.filling_count += 1,
            Event::None => {}
            _ => {}
        }
    }

    fn averaged(self, repetition_count: usize) -> AverageEventRow {
        debug_assert!(self.end_time > self.start_time);
        let repetitions = repetition_count as Float * EVENT_BIN_WIDTH as Float;
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

fn initial_event_row() -> AverageEventRow {
    AverageEventRow {
        time: 0.0,
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

fn temporary_result_paths(directory: &Path) -> Result<Vec<PathBuf>, String> {
    let entries = fs::read_dir(directory).map_err(|error| {
        format!(
            "failed to read temporary result directory {}: {error}",
            directory.display()
        )
    })?;
    let mut paths = Vec::new();

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
        if path.is_file() && name.starts_with("experiment_results_") && name.ends_with(".bin.gz") {
            paths.push(path);
        }
    }

    paths.sort();
    if paths.is_empty() {
        return Err(format!(
            "no Monte Carlo temporary result files were found in {}",
            directory.display()
        ));
    }
    Ok(paths)
}

/// Average the continuous fill trajectories in `tmp/` and write
/// `average_fill.csv` in the current experiment directory.
pub fn average_fill() -> Result<(), String> {
    average_fill_in(TEMPORARY_DIRECTORY, AVERAGE_FILL_FILE)
}

/// Average temporary trajectories from `temporary_directory` into
/// `output_file`.
pub fn average_fill_in(
    temporary_directory: impl AsRef<Path>,
    output_file: impl AsRef<Path>,
) -> Result<(), String> {
    let paths = temporary_result_paths(temporary_directory.as_ref())?;
    let rows = AverageFillRows::new(paths)?;
    write_continuous_values_csv(output_file, rows).map_err(|error| error.to_string())
}

/// Average event counts from `tmp/` into fixed-width bins and write
/// `average_event.csv` in the current experiment directory.
pub fn average_events() -> Result<(), String> {
    average_events_in(TEMPORARY_DIRECTORY, AVERAGE_EVENT_FILE)
}

/// Average temporary event trajectories from `temporary_directory` into
/// `output_file`.
pub fn average_events_in(
    temporary_directory: impl AsRef<Path>,
    output_file: impl AsRef<Path>,
) -> Result<(), String> {
    let paths = temporary_result_paths(temporary_directory.as_ref())?;
    let repetition_count = paths.len();
    let mut bins = Vec::<EventBin>::new();

    for path in paths {
        let batches =
            read_all_batches::<RecordedEvent>(&path).map_err(|error| error.to_string())?;
        for batch in batches {
            for record in batch.map_err(|error| error.to_string())? {
                validate_record(&path, &record)?;
                if record.time < 0.0 {
                    return Err(format!(
                        "{} contains a negative time {}, which cannot be placed in bins starting at zero",
                        path.display(),
                        record.time,
                    ));
                }

                let endpoint = (record.time / EVENT_BIN_WIDTH).ceil();
                if endpoint > usize::MAX as TimeFloat {
                    return Err(format!(
                        "{} contains a time too large to bin: {}",
                        path.display(),
                        record.time,
                    ));
                }
                let bins_through_record = endpoint as usize;
                while bins.len() < bins_through_record {
                    bins.push(EventBin::new(bins.len()));
                }

                if record.event != Event::None {
                    let bin_index = bins_through_record.saturating_sub(1);
                    if bins.is_empty() {
                        bins.push(EventBin::new(0));
                    }
                    bins[bin_index].record(record.event);
                }
            }
        }
    }

    let rows = std::iter::once(initial_event_row())
        .chain(bins.into_iter().map(|bin| bin.averaged(repetition_count)))
        .map(Ok::<_, std::convert::Infallible>);
    write_average_events_csv(output_file, rows).map_err(|error| error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::charge_transfer::Event;
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
        let output_file = directory.join("average_fill.csv");
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

        average_fill_in(&temporary_directory, &output_file).unwrap();
        let contents = fs::read_to_string(&output_file).unwrap();
        fs::remove_dir_all(directory).unwrap();

        let mut lines = contents.lines();
        assert_eq!(
            lines.next().unwrap(),
            "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9"
        );
        let rows = lines
            .map(|line| {
                line.split(',')
                    .map(|value| value.parse::<Float>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let expected = [
            [0.0, 100.0, 0.4, 0.2, 0.4, 0.24, 0.56],
            [1.0, 110.0, 0.5, 0.3, 0.5, 0.26, 0.74],
            [2.0, 120.0, 0.6, 0.2, 0.6, 0.44, 0.76],
            [3.0, 130.0, 0.3, 0.1, 0.3, 0.22, 0.38],
            [4.0, 140.0, 0.5, 0.1, 0.5, 0.42, 0.58],
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
        let output_file = directory.join("average_fill.csv");
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

        average_fill_in(&temporary_directory, &output_file).unwrap();
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
        let output_file = directory.join("average_event.csv");
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
                    Event::Filling {
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

        average_events_in(&temporary_directory, &output_file).unwrap();
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
        assert_eq!(rows[1][1], "15");
        assert_eq!(rows[1][3], "5");
        assert_eq!(rows[1][4], "5");
        assert_eq!(rows[1][9], "20");
        assert_eq!(rows[1][10], "5");
        assert_eq!(rows[1][11], "25");
        assert_eq!(rows[1][13], "5");
    }
}
