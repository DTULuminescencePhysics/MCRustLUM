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
use std::error::Error;
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

/// An error produced while writing or reading temporary experiment output.
#[derive(Debug)]
pub enum OutputError {
    /// The temporary output file could not be opened.
    Open {
        /// File that could not be opened.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
    /// A batch could not be encoded or written.
    Write {
        /// File to which serialization failed.
        path: PathBuf,
        /// Bincode serialization or output error.
        source: Box<bincode::ErrorKind>,
    },
    /// The final bytes of a gzip member could not be written.
    Finish {
        /// File whose gzip member could not be finalized.
        path: PathBuf,
        /// Underlying compression or filesystem error.
        source: std::io::Error,
    },
    /// Compressed data could not be read.
    Read {
        /// Compressed file that could not be read.
        path: PathBuf,
        /// Underlying decompression or filesystem error.
        source: std::io::Error,
    },
    /// A batch was incomplete or did not match the requested record type.
    Decode {
        /// File containing the malformed or type-incompatible batch.
        path: PathBuf,
        /// Bincode deserialization error.
        source: Box<bincode::ErrorKind>,
    },
}

impl fmt::Display for OutputError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Open { path, source } => {
                write!(formatter, "failed to open {}: {source}", path.display())
            }
            Self::Write { path, source } => {
                write!(
                    formatter,
                    "failed to write batch to {}: {source}",
                    path.display()
                )
            }
            Self::Finish { path, source } => write!(
                formatter,
                "failed to finish compressed batch in {}: {source}",
                path.display(),
            ),
            Self::Read { path, source } => {
                write!(formatter, "failed to read {}: {source}", path.display())
            }
            Self::Decode { path, source } => {
                write!(
                    formatter,
                    "failed to decode batch from {}: {source}",
                    path.display()
                )
            }
        }
    }
}

impl Error for OutputError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Open { source, .. } | Self::Finish { source, .. } | Self::Read { source, .. } => {
                Some(source)
            }
            Self::Write { source, .. } | Self::Decode { source, .. } => Some(source.as_ref()),
        }
    }
}

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
/// # fn example() -> Result<(), io::outputs::OutputError> {
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
#[derive(Debug, Clone, Copy, PartialEq)]
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
}

/// One row in the averaged event-count output.
///
/// Every count is normalized by both repetition count and bin width and is
/// therefore an observed event frequency in s⁻¹, not a microscopic rate
/// equation evaluated at the row time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AverageEventRow {
    /// Right-hand edge of the event bin, in seconds.
    pub time: TimeFloat,
    /// Ground-state localised recombinations per second per repetition.
    pub localised_recombination_ground_count: Float,
    /// Excited-state localised recombinations per second per repetition.
    pub localised_recombination_excited_count: Float,
    /// Ground-state delocalised recombinations per second per repetition.
    pub delocalised_recombination_ground_count: Float,
    /// Excited-state delocalised recombinations per second per repetition.
    pub delocalised_recombination_excited_count: Float,
    /// Ground-state localised retrapping events per second per repetition.
    pub localised_retrapping_ground_count: Float,
    /// Excited-state localised retrapping events per second per repetition.
    pub localised_retrapping_excited_count: Float,
    /// Ground-state delocalised retrapping events per second per repetition.
    pub delocalised_retrapping_ground_count: Float,
    /// Excited-state delocalised retrapping events per second per repetition.
    pub delocalised_retrapping_excited_count: Float,
    /// All events originating from ground states per second per repetition.
    pub ground_count: Float,
    /// All events originating from excited states per second per repetition.
    pub excited_count: Float,
    /// All recombination events per second per repetition.
    pub recombination_count: Float,
    /// All retrapping events per second per repetition.
    pub retrapping_count: Float,
    /// All irradiation-driven filling events per second per repetition.
    pub filling_count: Float,
}

/// An error produced while writing consolidated CSV output.
#[derive(Debug)]
pub enum CsvOutputError {
    /// The destination CSV file could not be created.
    Create {
        /// Destination path.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
    /// The input iterator could not produce a valid output row.
    SourceData {
        /// Destination path for which rows were being generated.
        path: PathBuf,
        /// Display form of the upstream row-generation error.
        message: String,
    },
    /// A header, row, or buffered tail could not be written.
    Write {
        /// Destination path.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
}

impl fmt::Display for CsvOutputError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Create { path, source } => {
                write!(formatter, "failed to create {}: {source}", path.display())
            }
            Self::SourceData { path, message } => write!(
                formatter,
                "failed to produce a row for {}: {message}",
                path.display()
            ),
            Self::Write { path, source } => {
                write!(formatter, "failed to write {}: {source}", path.display())
            }
        }
    }
}

impl Error for CsvOutputError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Create { source, .. } | Self::Write { source, .. } => Some(source),
            Self::SourceData { .. } => None,
        }
    }
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
        "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9"
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
            "{},{},{},{},{},{},{}",
            row.time,
            row.temperature,
            row.fill,
            row.fill_standard_deviation,
            row.fill_median,
            row.fill_quantile_0_1,
            row.fill_quantile_0_9,
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

/// Write normalized event-frequency rows to a CSV file.
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
        "time,localised_recombination_ground,localised_recombination_excited,delocalised_recombination_ground,delocalised_recombination_excited,localised_retrapping_ground,localised_retrapping_excited,delocalised_retrapping_ground,delocalised_retrapping_excited,ground,excited,recombination, retrapping, filling_count"
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
            },
            ContinuousValueRow {
                time: 1.0,
                temperature: 283.15,
                fill: 0.5,
                fill_standard_deviation: 0.1,
                fill_median: 0.5,
                fill_quantile_0_1: 0.42,
                fill_quantile_0_9: 0.58,
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
                "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9\n",
                "0,273.15,0.25,0.05,0.25,0.21,0.29\n",
                "1,283.15,0.5,0.1,0.5,0.42,0.58\n",
            )
        );
    }
}
