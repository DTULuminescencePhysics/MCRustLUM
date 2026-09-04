// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Consolidation of repeated Monte Carlo trajectories.

use common::charge_transfer::RecordedEvent;
use common::numeric::{Float, TimeFloat};
use io::outputs::{BatchReader, ContinuousValueRow, read_all_batches, write_continuous_values_csv};
use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};
use std::vec::IntoIter;

const TEMPORARY_DIRECTORY: &str = "tmp";
const AVERAGE_FILL_FILE: &str = "average_fill.csv";
const AVERAGE_EVENT_FILE: &str = "average_fill.csv";

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
    fill_sum: Float,
    first: Option<ContinuousValueRow>,
    finished: bool,
}

impl AverageFillRows {
    fn new(paths: Vec<PathBuf>) -> Result<Self, String> {
        let mut cursors = Vec::with_capacity(paths.len());
        let mut initial_time = None;
        let mut initial_temperature_sum = 0.0;
        let mut fill_sum = 0.0;

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
            fill_sum += first.fill;
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
        let first = ContinuousValueRow {
            time: initial_time.expect("at least one path is required"),
            temperature: initial_temperature_sum / count,
            fill: fill_sum / count,
        };

        Ok(Self {
            cursors,
            direction,
            fill_sum,
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
        let mut fill_sum = 0.0;
        for cursor in &mut self.cursors {
            match Self::advance_cursor_at(cursor, time, self.direction) {
                Ok((cursor_temperature_sum, cursor_temperature_count)) => {
                    fill_sum += cursor.fill;
                    temperature_sum += cursor_temperature_sum;
                    temperature_count += cursor_temperature_count;
                }
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            }
        }
        self.fill_sum = fill_sum;

        if temperature_count == 0 {
            self.finished = true;
            return Some(Err(format!("no temperature was recorded at time {time}")));
        }

        Some(Ok(ContinuousValueRow {
            time,
            temperature: temperature_sum / temperature_count as Float,
            fill: self.fill_sum / self.cursors.len() as Float,
        }))
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
    
    filling_count: usize

}

struct AverageEventRows {
    cursors: Vec<TrajectoryCursor>,
    direction: Direction,
    fill_sum: Float,
    first: Option<ContinuousValueRow>,
    finished: bool,

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

#[cfg(test)]
mod tests {
    use super::*;
    use common::charge_transfer::Event;
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

        assert_eq!(
            contents,
            concat!(
                "time,temperature,fill\n",
                "0,100,0.4\n",
                "1,110,0.5\n",
                "2,120,0.6000000000000001\n",
                "3,130,0.30000000000000004\n",
                "4,140,0.5\n",
            )
        );
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
}
