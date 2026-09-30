// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Streaming output for chronology profiles.
//!
//! The uncompressed bincode stream contains, in order:
//!
//! 1. one [`TimeUnit`];
//! 2. one [`TemperatureUnit`];
//! 3. zero or more `Vec<T>` batches.
//!
//! The header and every appended batch are separate gzip members. A
//! [`MultiGzDecoder`] exposes them as one continuous bincode stream when read.

use crate::errors::{CsvOutputError, OutputError};
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bincode::Options;
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde::de::DeserializeOwned;
use serde::Serialize;

/// Return the stable bincode configuration shared by the writer and reader.
fn bincode_options() -> impl Options {
    bincode::DefaultOptions::new().with_fixint_encoding()
}

/// Create an empty chronology output file, replacing any previous file.
///
/// Call [`write_temporary_chronology_units`] next, before appending batches.
pub fn create_temporary_chronology_experiment_file(
    path: impl AsRef<Path>,
) -> Result<(), OutputError> {
    let path = path.as_ref();
    File::create(path).map_err(|source| OutputError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

/// Write the two-entry unit header and discard any old contents of the file.
///
/// Truncating here guarantees that the units are always the first two values
/// in the bincode stream, even if a previous run left a file at `path`.
pub fn write_temporary_chronology_units(
    path: impl AsRef<Path>,
    time_unit: TimeUnit,
    temp_unit: TemperatureUnit,
) -> Result<(), OutputError> {
    let path = path.as_ref();
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(path)
        .map_err(|source| OutputError::Open {
            path: path.to_path_buf(),
            source,
        })?;

    let writer = BufWriter::new(file);
    let mut encoder = GzEncoder::new(writer, Compression::default());

    bincode_options()
        .serialize_into(&mut encoder, &time_unit)
        .map_err(|source| OutputError::Write {
            path: path.to_path_buf(),
            source,
        })?;
    bincode_options()
        .serialize_into(&mut encoder, &temp_unit)
        .map_err(|source| OutputError::Write {
            path: path.to_path_buf(),
            source,
        })?;

    let mut writer = encoder.finish().map_err(|source| OutputError::Finish {
        path: path.to_path_buf(),
        source,
    })?;
    writer.flush().map_err(|source| OutputError::Finish {
        path: path.to_path_buf(),
        source,
    })
}

/// Append one vector of profiles as a new compressed batch.
///
/// `T` is generic so this crate does not need to depend on the `chronology`
/// crate (which already depends on this crate). In normal use, `T` is
/// `chronology::profiles::TimeTempProfile`.
pub fn append_chronology_profile_to_experiment_file<T: Serialize>(
    path: impl AsRef<Path>,
    batch: &[T],
) -> Result<(), OutputError> {
    if batch.is_empty() {
        return Ok(());
    }

    let path = path.as_ref();
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .map_err(|source| OutputError::Open {
            path: path.to_path_buf(),
            source,
        })?;
    let writer = BufWriter::new(file);
    let mut encoder = GzEncoder::new(writer, Compression::default());

    bincode_options()
        .serialize_into(&mut encoder, batch)
        .map_err(|source| OutputError::Write {
            path: path.to_path_buf(),
            source,
        })?;

    let mut writer = encoder.finish().map_err(|source| OutputError::Finish {
        path: path.to_path_buf(),
        source,
    })?;
    writer.flush().map_err(|source| OutputError::Finish {
        path: path.to_path_buf(),
        source,
    })
}

/// Buffered decoder that treats concatenated gzip members as one byte stream.
type CompressedReader = BufReader<MultiGzDecoder<BufReader<File>>>;

/// The unit header and a streaming iterator over the stored profile batches.
pub struct ChronologyBatchReader<T> {
    /// Unit used by every stored time value.
    pub time_unit: TimeUnit,
    /// Unit used by every stored temperature value.
    pub temp_unit: TemperatureUnit,
    path: PathBuf,
    reader: CompressedReader,
    finished: bool,
    record: PhantomData<T>,
}

impl<T: DeserializeOwned> Iterator for ChronologyBatchReader<T> {
    type Item = Result<Vec<T>, OutputError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        match self.reader.fill_buf() {
            Ok(buffer) if buffer.is_empty() => {
                self.finished = true;
                None
            }
            Ok(_) => match bincode_options().deserialize_from(&mut self.reader) {
                Ok(batch) => Some(Ok(batch)),
                Err(source) => {
                    self.finished = true;
                    Some(Err(OutputError::Decode {
                        path: self.path.clone(),
                        source,
                    }))
                }
            },
            Err(source) => {
                self.finished = true;
                Some(Err(OutputError::Read {
                    path: self.path.clone(),
                    source,
                }))
            }
        }
    }
}

/// Open a chronology file, read its unit header, and stream its batches.
pub fn read_all_chrono_batches<T: DeserializeOwned>(
    path: impl AsRef<Path>,
) -> Result<ChronologyBatchReader<T>, OutputError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|source| OutputError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    let compressed = MultiGzDecoder::new(BufReader::new(file));
    let mut reader = BufReader::new(compressed);

    let time_unit = bincode_options()
        .deserialize_from(&mut reader)
        .map_err(|source| OutputError::Decode {
            path: path.to_path_buf(),
            source,
        })?;
    let temp_unit = bincode_options()
        .deserialize_from(&mut reader)
        .map_err(|source| OutputError::Decode {
            path: path.to_path_buf(),
            source,
        })?;

    Ok(ChronologyBatchReader {
        time_unit,
        temp_unit,
        path: path.to_path_buf(),
        reader,
        finished: false,
        record: PhantomData,
    })
}

pub struct ProfileGrid {
    pub accepted_profiles: u64,

    pub time_bounds: (TimeFloat, TimeFloat),
    pub temperature_bounds: (Float, Float),

    pub time_unit: TimeUnit,
    pub temperature_unit: TemperatureUnit,

    pub time_cells: usize,
    pub temperature_cells: usize,

    /// Row-major: `weights[temp_index * time_cells + time_index]`.
    ///
    /// Both indices are stored internally in ascending physical order.
    pub weights: Vec<u64>,
}

pub fn write_profile_grid_csv(
    path: impl AsRef<Path>,
    grid: &ProfileGrid,
) -> Result<(), CsvOutputError> {
    let path = path.as_ref();
    let expected_cells = grid
        .time_cells
        .checked_mul(grid.temperature_cells)
        .filter(|&cells| cells > 0)
        .ok_or_else(|| CsvOutputError::SourceData {
            path: path.to_path_buf(),
            message: "grid dimensions must be non-zero and must not overflow".to_string(),
        })?;
    if grid.weights.len() != expected_cells {
        return Err(CsvOutputError::SourceData {
            path: path.to_path_buf(),
            message: format!(
                "grid has {} weights but its dimensions require {expected_cells}",
                grid.weights.len()
            ),
        });
    }

    let file = File::create(path).map_err(|source| CsvOutputError::Create {
        path: path.to_path_buf(),
        source,
    })?;

    let mut writer = BufWriter::new(file);

    // First line: accepted-profile count.
    writeln!(writer, "{}", grid.accepted_profiles).map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })?;

    let reverse_time = grid.time_unit.is_ka_or_ma();

    let (time_start, time_end) = if reverse_time {
        (grid.time_bounds.1, grid.time_bounds.0)
    } else {
        (grid.time_bounds.0, grid.time_bounds.1)
    };

    // Second line: start time, end time, unit.
    writeln!(writer, "{time_start},{time_end},{:?}", grid.time_unit).map_err(|source| {
        CsvOutputError::Write {
            path: path.to_path_buf(),
            source,
        }
    })?;

    // Third line: minimum temperature, maximum temperature, unit.
    writeln!(
        writer,
        "{},{},{:?}",
        grid.temperature_bounds.0, grid.temperature_bounds.1, grid.temperature_unit
    )
    .map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })?;

    // Highest-temperature row first, as with a conventional y-axis.
    for temperature_index in (0..grid.temperature_cells).rev() {
        for output_column in 0..grid.time_cells {
            if output_column != 0 {
                write!(writer, ",").map_err(|source| CsvOutputError::Write {
                    path: path.to_path_buf(),
                    source,
                })?;
            }

            // Geological time runs t -> 0, so reverse its columns.
            let time_index = if reverse_time {
                grid.time_cells - 1 - output_column
            } else {
                output_column
            };

            let index = temperature_index * grid.time_cells + time_index;
            write!(writer, "{}", grid.weights[index]).map_err(|source| CsvOutputError::Write {
                path: path.to_path_buf(),
                source,
            })?;
        }

        writeln!(writer).map_err(|source| CsvOutputError::Write {
            path: path.to_path_buf(),
            source,
        })?;
    }

    writer.flush().map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::{Deserialize, Serialize};
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[derive(Debug, PartialEq, Serialize, Deserialize)]
    struct TestProfile {
        times: Vec<f64>,
        temps: Vec<f64>,
    }

    fn temporary_output_path() -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_chronology_{}_{}.bin.gz",
            std::process::id(),
            unique,
        ))
    }

    #[test]
    fn writes_units_once_and_streams_variable_length_batches() {
        let path = temporary_output_path();
        let first = vec![TestProfile {
            times: vec![0.0, 1.0],
            temps: vec![280.0, 281.0],
        }];
        let second = vec![
            TestProfile {
                times: vec![2.0],
                temps: vec![282.0],
            },
            TestProfile {
                times: vec![3.0, 4.0, 5.0],
                temps: vec![283.0, 284.0, 285.0],
            },
        ];

        write_temporary_chronology_units(&path, TimeUnit::KAnnum, TemperatureUnit::Celsius)
            .unwrap();
        append_chronology_profile_to_experiment_file(&path, &first).unwrap();
        append_chronology_profile_to_experiment_file(&path, &second).unwrap();

        let reader = read_all_chrono_batches::<TestProfile>(&path).unwrap();
        assert_eq!(reader.time_unit, TimeUnit::KAnnum);
        assert_eq!(reader.temp_unit, TemperatureUnit::Celsius);
        let batches = reader.collect::<Result<Vec<_>, _>>().unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(batches, vec![first, second]);
    }

    #[test]
    fn writes_grid_in_display_axis_order() {
        let path = temporary_output_path();
        let grid = ProfileGrid {
            accepted_profiles: 7,
            time_bounds: (0.0, 10.0),
            temperature_bounds: (20.0, 40.0),
            time_unit: TimeUnit::KAnnum,
            temperature_unit: TemperatureUnit::Celsius,
            time_cells: 2,
            temperature_cells: 2,
            // Ascending temperature rows and ascending time columns.
            weights: vec![1, 2, 3, 4],
        };

        write_profile_grid_csv(&path, &grid).unwrap();
        let output = fs::read_to_string(&path).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(output, "7\n10,0,KAnnum\n20,40,Celsius\n4,3\n2,1\n");
    }

    #[test]
    fn rejects_inconsistent_grid_dimensions() {
        let path = temporary_output_path();
        let grid = ProfileGrid {
            accepted_profiles: 0,
            time_bounds: (0.0, 1.0),
            temperature_bounds: (0.0, 1.0),
            time_unit: TimeUnit::Second,
            temperature_unit: TemperatureUnit::Kelvin,
            time_cells: 2,
            temperature_cells: 2,
            weights: vec![0; 3],
        };

        let error = write_profile_grid_csv(&path, &grid).unwrap_err();

        assert!(matches!(error, CsvOutputError::SourceData { .. }));
        assert!(!path.exists());
    }
}
