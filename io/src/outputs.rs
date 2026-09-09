// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Streaming temporary output for Monte Carlo experiments.
//!
//! Each call to [`crate::outputs::append_monte_carlo_experiment_batch_to_file`]
//! adds one bincode-encoded batch as a new gzip member.
//! [`crate::outputs::read_all_batches`] decodes
//! those members as a stream, keeping only the current batch in memory.

use common::numeric::{Float, TimeFloat};
use crate::errors::{OutputError, CsvOutputError};
use std::fmt;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bincode::Options;
use flate2::Compression;
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use serde::Serialize;
use serde::de::DeserializeOwned;


/// Return the stable bincode configuration shared by the writer and reader.
fn bincode_options() -> impl Options {
    bincode::DefaultOptions::new().with_fixint_encoding()
}

/// Create an empty temporary experiment file, replacing a previous file at
/// the same path.
///
/// Call this once before appending the first batch of a new experiment. This
/// prevents a repeated simulation run from accidentally retaining old batches.
pub fn create_monte_carlo_experiment_file(path: impl AsRef<Path>) -> Result<(), OutputError> {
    let path = path.as_ref();
    File::create(path).map_err(|source| OutputError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

/// Append one batch of records to a gzip-compressed temporary file.
///
/// The records are serialized directly into the compressor. An empty batch is
/// ignored. The caller can therefore keep a bounded `Vec`, call this function
/// with a slice, and then clear and reuse that allocation.
pub fn append_monte_carlo_experiment_batch_to_file<T: Serialize>(
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
    })?;

    Ok(())
}

/// Buffered decoder that treats concatenated gzip members as one byte stream.
type CompressedReader = BufReader<MultiGzDecoder<BufReader<File>>>;

/// Iterator over the batches stored in one temporary experiment file.
///
/// A batch is dropped before the next one needs to be decoded, provided the
/// caller does not retain it. After the first read or decode error the iterator
/// is exhausted.
pub struct BatchReader<T> {
    /// Source path retained for contextual read and decode errors.
    path: PathBuf,
    /// Streaming decoder positioned at the next bincode batch.
    reader: CompressedReader,
    /// Whether end-of-file or an unrecoverable error has been observed.
    finished: bool,
    /// Associates the iterator with the deserialized record type.
    record: PhantomData<T>,
}

impl<T: DeserializeOwned> Iterator for BatchReader<T> {
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

/// Open a temporary experiment file and stream all of its batches.
///
/// # Example
///
/// ```no_run
/// # use serde::Deserialize;
/// # #[derive(Deserialize)]
/// # struct Event;
/// # fn process(_: Vec<Event>) {}
/// # fn example() -> Result<(), io::errors::OutputError> {
/// for batch in io::outputs::read_all_batches::<Event>("experiment.bin.gz")? {
///     process(batch?);
/// }
/// # Ok(())
/// # }
/// ```
pub fn read_all_batches<T: DeserializeOwned>(
    path: impl AsRef<Path>,
) -> Result<BatchReader<T>, OutputError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|source| OutputError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    let compressed = MultiGzDecoder::new(BufReader::new(file));

    Ok(BatchReader {
        path: path.to_path_buf(),
        reader: BufReader::new(compressed),
        finished: false,
        record: PhantomData,
    })
}

/// One row in the consolidated continuous-value output.
///
/// Fill statistics are calculated across independent repetitions at the same
/// profile time. Temperature is the mean of records that occur at that time.
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
pub struct ContinuousValueRow {
    /// Absolute profile time in seconds.
    pub time: TimeFloat,
    /// Mean temperature of contributing records, in kelvin.
    pub temperature: Float,
    /// Arithmetic mean trap filling fraction.
    pub fill: Float,
    /// Population standard deviation of the filling fraction.
    pub fill_standard_deviation: Float,
    /// Median filling fraction, using linearly interpolated quantiles.
    pub fill_median: Float,
    /// Linearly interpolated 10th percentile of the filling fraction.
    pub fill_quantile_0_1: Float,
    /// Linearly interpolated 90th percentile of the filling fraction.
    pub fill_quantile_0_9: Float,
    /// Linearly interpolated lower quartile of the filling fraction.
    #[serde(default = "missing_quantile")]
    pub fill_quantile_0_25: Float,
    /// Linearly interpolated upper quartile of the filling fraction.
    #[serde(default = "missing_quantile")]
    pub fill_quantile_0_75: Float,
}

/// Mark quartiles as unavailable when reading CSV files made by older releases.
fn missing_quantile() -> Float {
    Float::NAN
}

/// One row in the averaged event-count output.
///
/// Every count is normalized by the repetition count. Divide a value by the
/// duration between this row's right edge and the preceding edge to obtain an
/// observed event frequency in s⁻¹.
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
pub struct AverageEventRow {
    /// Right-hand edge of the event bin, in seconds.
    pub time: TimeFloat,
    /// Ground-state localised recombinations per repetition in this bin.
    #[serde(rename = "localised_recombination_ground")]
    pub localised_recombination_ground_count: Float,
    /// Excited-state localised recombinations per repetition in this bin.
    #[serde(rename = "localised_recombination_excited")]
    pub localised_recombination_excited_count: Float,
    /// Ground-state delocalised recombinations per repetition in this bin.
    #[serde(rename = "delocalised_recombination_ground")]
    pub delocalised_recombination_ground_count: Float,
    /// Excited-state delocalised recombinations per repetition in this bin.
    #[serde(rename = "delocalised_recombination_excited")]
    pub delocalised_recombination_excited_count: Float,
    /// Ground-state localised retrapping events per repetition in this bin.
    #[serde(rename = "localised_retrapping_ground")]
    pub localised_retrapping_ground_count: Float,
    /// Excited-state localised retrapping events per repetition in this bin.
    #[serde(rename = "localised_retrapping_excited")]
    pub localised_retrapping_excited_count: Float,
    /// Ground-state delocalised retrapping events per repetition in this bin.
    #[serde(rename = "delocalised_retrapping_ground")]
    pub delocalised_retrapping_ground_count: Float,
    /// Excited-state delocalised retrapping events per repetition in this bin.
    #[serde(rename = "delocalised_retrapping_excited")]
    pub delocalised_retrapping_excited_count: Float,
    /// All events originating from ground states per repetition in this bin.
    #[serde(rename = "ground")]
    pub ground_count: Float,
    /// All events originating from excited states per repetition in this bin.
    #[serde(rename = "excited")]
    pub excited_count: Float,
    /// All recombination events per repetition in this bin.
    #[serde(rename = "recombination")]
    pub recombination_count: Float,
    /// All retrapping events per repetition in this bin.
    #[serde(rename = "retrapping")]
    pub retrapping_count: Float,
    /// All irradiation-driven filling events per repetition in this bin.
    pub filling_count: Float,
}



/// Write consolidated time, temperature, and fill rows to a CSV file.
///
/// Rows are consumed one at a time, allowing the caller to consolidate large
/// temporary files without first collecting the final output in memory.
pub fn write_continuous_values_csv<I, E>(
    path: impl AsRef<Path>,
    rows: I,
) -> Result<(), CsvOutputError>
where
    I: IntoIterator<Item = Result<ContinuousValueRow, E>>,
    E: fmt::Display,
{
    let path = path.as_ref();
    let file = File::create(path).map_err(|source| CsvOutputError::Create {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9,fill_quantile_0_25,fill_quantile_0_75"
    )
    .map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })?;

    for row in rows {
        let row = row.map_err(|error| CsvOutputError::SourceData {
            path: path.to_path_buf(),
            message: error.to_string(),
        })?;
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{}",
            row.time,
            row.temperature,
            row.fill,
            row.fill_standard_deviation,
            row.fill_median,
            row.fill_quantile_0_1,
            row.fill_quantile_0_9,
            row.fill_quantile_0_25,
            row.fill_quantile_0_75,
        )
        .map_err(|source| CsvOutputError::Write {
            path: path.to_path_buf(),
            source,
        })?;
    }

    writer.flush().map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })
}

/// Write per-bin event counts averaged across repetitions to a CSV file.
///
/// Rows are streamed and preserve the category breakdown in
/// [`AverageEventRow`]. The function writes column names before requesting the
/// first row, so an upstream error may leave a header-only partial file.
pub fn write_average_events_csv<I, E>(path: impl AsRef<Path>, rows: I) -> Result<(), CsvOutputError>
where
    I: IntoIterator<Item = Result<AverageEventRow, E>>,
    E: fmt::Display,
{
    let path = path.as_ref();
    let file = File::create(path).map_err(|source| CsvOutputError::Create {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "time,localised_recombination_ground,localised_recombination_excited,delocalised_recombination_ground,delocalised_recombination_excited,localised_retrapping_ground,localised_retrapping_excited,delocalised_retrapping_ground,delocalised_retrapping_excited,ground,excited,recombination,retrapping,filling_count"
    )
    .map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })?;

    for row in rows {
        let row = row.map_err(|error| CsvOutputError::SourceData {
            path: path.to_path_buf(),
            message: error.to_string(),
        })?;
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            row.time,
            row.localised_recombination_ground_count,
            row.localised_recombination_excited_count,
            row.delocalised_recombination_ground_count,
            row.delocalised_recombination_excited_count,
            row.localised_retrapping_ground_count,
            row.localised_retrapping_excited_count,
            row.delocalised_retrapping_ground_count,
            row.delocalised_retrapping_excited_count,
            row.ground_count,
            row.excited_count,
            row.recombination_count,
            row.retrapping_count,
            row.filling_count,
        )
        .map_err(|source| CsvOutputError::Write {
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
    use std::convert::Infallible;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};
    #[derive(Debug, PartialEq, Serialize, Deserialize)]
    struct TestRecord {
        id: u16,
        value: f64,
    }

    fn temporary_output_path(label: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_{label}_{}_{}.bin.gz",
            std::process::id(),
            unique,
        ))
    }

    #[test]
    fn appends_and_streams_multiple_batches() {
        let path = temporary_output_path("batches");
        let first = vec![
            TestRecord { id: 1, value: 1.5 },
            TestRecord { id: 2, value: 2.5 },
        ];
        let second = vec![TestRecord { id: 3, value: 3.5 }];

        append_monte_carlo_experiment_batch_to_file(&path, &first).unwrap();
        append_monte_carlo_experiment_batch_to_file(&path, &second).unwrap();

        let batches = read_all_batches::<TestRecord>(&path)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        fs::remove_file(&path).unwrap();

        assert_eq!(batches, vec![first, second]);
    }

    #[test]
    fn an_empty_batch_does_not_create_a_file() {
        let path = temporary_output_path("empty");

        append_monte_carlo_experiment_batch_to_file::<TestRecord>(&path, &[]).unwrap();

        assert!(!path.exists());
    }

    #[test]
    fn creating_an_experiment_file_discards_old_batches() {
        let path = temporary_output_path("truncate");
        let old = vec![TestRecord { id: 1, value: 1.5 }];
        let replacement = vec![TestRecord { id: 2, value: 2.5 }];

        append_monte_carlo_experiment_batch_to_file(&path, &old).unwrap();
        create_monte_carlo_experiment_file(&path).unwrap();
        append_monte_carlo_experiment_batch_to_file(&path, &replacement).unwrap();

        let batches = read_all_batches::<TestRecord>(&path)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        fs::remove_file(&path).unwrap();

        assert_eq!(batches, vec![replacement]);
    }

    #[test]
    fn missing_input_reports_the_path() {
        let path = temporary_output_path("missing");
        let error = match read_all_batches::<TestRecord>(&path) {
            Ok(_) => panic!("missing temporary output should fail"),
            Err(error) => error,
        };

        assert!(error.to_string().contains(path.to_string_lossy().as_ref()));
        assert!(matches!(error, OutputError::Open { .. }));
    }

    fn temporary_output_path_2() -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_continuous_{}_{}.csv",
            std::process::id(),
            unique,
        ))
    }

    #[test]
    fn writes_header_and_continuous_rows() {
        let path = temporary_output_path_2();
        let rows = [
            ContinuousValueRow {
                time: 0.0,
                temperature: 273.15,
                fill: 0.25,
                fill_standard_deviation: 0.05,
                fill_median: 0.25,
                fill_quantile_0_1: 0.21,
                fill_quantile_0_9: 0.29,
                fill_quantile_0_25: 0.225,
                fill_quantile_0_75: 0.275,
            },
            ContinuousValueRow {
                time: 1.0,
                temperature: 283.15,
                fill: 0.5,
                fill_standard_deviation: 0.1,
                fill_median: 0.5,
                fill_quantile_0_1: 0.42,
                fill_quantile_0_9: 0.58,
                fill_quantile_0_25: 0.45,
                fill_quantile_0_75: 0.55,
            },
        ]
        .into_iter()
        .map(Ok::<_, Infallible>);

        write_continuous_values_csv(&path, rows).unwrap();
        let contents = fs::read_to_string(&path).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(
            contents,
            concat!(
                "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9,fill_quantile_0_25,fill_quantile_0_75\n",
                "0,273.15,0.25,0.05,0.25,0.21,0.29,0.225,0.275\n",
                "1,283.15,0.5,0.1,0.5,0.42,0.58,0.45,0.55\n",
            )
        );
    }
}
