// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Streaming temporary output for Monte Carlo experiments.
//!
//! Each call to [`append_monte_carlo_experiment_batch_to_file`] adds one
//! bincode-encoded batch as a new gzip member. [`read_all_batches`] decodes
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
        path: PathBuf,
        source: std::io::Error,
    },
    /// A batch could not be encoded or written.
    Write {
        path: PathBuf,
        source: Box<bincode::ErrorKind>,
    },
    /// The final bytes of a gzip member could not be written.
    Finish {
        path: PathBuf,
        source: std::io::Error,
    },
    /// Compressed data could not be read.
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    /// A batch was incomplete or did not match the requested record type.
    Decode {
        path: PathBuf,
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

type CompressedReader = BufReader<MultiGzDecoder<BufReader<File>>>;

/// Iterator over the batches stored in one temporary experiment file.
///
/// A batch is dropped before the next one needs to be decoded, provided the
/// caller does not retain it. After the first read or decode error the iterator
/// is exhausted.
pub struct BatchReader<T> {
    path: PathBuf,
    reader: CompressedReader,
    finished: bool,
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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContinuousValueRow {
    pub time: TimeFloat,
    pub temperature: Float,
    pub fill: Float,
}

/// An error produced while writing consolidated CSV output.
#[derive(Debug)]
pub enum CsvOutputError {
    Create {
        path: PathBuf,
        source: std::io::Error,
    },
    SourceData {
        path: PathBuf,
        message: String,
    },
    Write {
        path: PathBuf,
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

    writeln!(writer, "time,temperature,fill").map_err(|source| CsvOutputError::Write {
        path: path.to_path_buf(),
        source,
    })?;

    for row in rows {
        let row = row.map_err(|error| CsvOutputError::SourceData {
            path: path.to_path_buf(),
            message: error.to_string(),
        })?;
        writeln!(writer, "{},{},{}", row.time, row.temperature, row.fill).map_err(|source| {
            CsvOutputError::Write {
                path: path.to_path_buf(),
                source,
            }
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
    use std::convert::Infallible;
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
            },
            ContinuousValueRow {
                time: 1.0,
                temperature: 283.15,
                fill: 0.5,
            },
        ]
        .into_iter()
        .map(Ok::<_, Infallible>);

        write_continuous_values_csv(&path, rows).unwrap();
        let contents = fs::read_to_string(&path).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(
            contents,
            "time,temperature,fill\n0,273.15,0.25\n1,283.15,0.5\n"
        );
    }

}
