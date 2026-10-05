// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Error catching for the filesystem and input and outputs

use std::error::Error;
use std::fmt;
use std::path::PathBuf;

/// An error produced while preparing an experiment directory.
#[derive(Debug)]
pub enum FilesystemError {
    /// A supplied folder name was empty, nested, or otherwise unsafe.
    InvalidFolderName {
        /// Rejected path supplied as the experiment name.
        name: PathBuf,
    },
    /// A filesystem operation failed.
    Operation {
        /// Human-readable filesystem operation being attempted.
        action: &'static str,
        /// Path on which the operation failed.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
    /// Every representable automatic experiment number was already occupied.
    ExperimentNumberExhausted,
}

impl fmt::Display for FilesystemError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidFolderName { name } => write!(
                formatter,
                "experiment folder name must be one normal path component, got {:?}",
                name
            ),
            Self::Operation {
                action,
                path,
                source,
            } => write!(formatter, "failed to {action} {}: {source}", path.display()),
            Self::ExperimentNumberExhausted => {
                formatter.write_str("could not find an available automatic experiment number")
            }
        }
    }
}

impl Error for FilesystemError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Operation { source, .. } => Some(source),
            Self::InvalidFolderName { .. } | Self::ExperimentNumberExhausted => None,
        }
    }
}

/// An error produced while reading or parsing a simulation input file.
#[derive(Debug)]
pub enum InputError {
    /// The input file could not be opened or read.
    Read {
        /// Path that could not be read.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
    /// File contents were not valid TOML for [`crate::inputs::SimulationInputs`]
    /// or for [`crate::chrono_inputs::ChronologyInputs`].
    Parse {
        /// Path containing the invalid TOML.
        path: PathBuf,
        /// Underlying TOML deserialization error.
        source: toml::de::Error,
    },
}
impl fmt::Display for InputError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Read { path, source } => {
                write!(formatter, "failed to read {}: {source}", path.display())
            }
            Self::Parse { path, source } => {
                write!(formatter, "failed to parse {}: {source}", path.display())
            }
        }
    }
}
impl Error for InputError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Read { source, .. } => Some(source),
            Self::Parse { source, .. } => Some(source),
        }
    }
}
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
    Invalid {
        /// File attempting to be written to.
        path: PathBuf,
        /// Bincode serialization or output error.
        source: Box<bincode::ErrorKind>,
        /// Human-readable filesystem operation being attempted.
        action: &'static str,
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
            Self::Invalid {
                path,
                source,
                action,
            } => {
                write!(
                    formatter,
                    "Error message: {action}. So failed to write batch to {}: {source}",
                    path.display(),
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
            Self::Write { source, .. }
            | Self::Decode { source, .. }
            | Self::Invalid { source, .. } => Some(source.as_ref()),
        }
    }
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

/// Error returned while loading, validating, rebinning, or plotting results.
#[derive(Debug)]
pub enum PlotError {
    Setup {
        source: String,
        message: String,
    },
    /// A CSV file could not be opened or parsed.
    Csv {
        /// Input path associated with the error.
        path: PathBuf,
        /// Error reported by the CSV reader.
        source: csv::Error,
    },
    /// Result data was empty or physically inconsistent.
    InvalidData {
        /// Input path containing the invalid data.
        path: PathBuf,
        /// Explanation of the violated requirement.
        message: String,
    },
    /// An output directory could not be created.
    CreateDirectory {
        /// Directory that could not be created.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
    /// Plotters could not draw or save an image.
    Draw {
        /// Destination image path.
        path: PathBuf,
        /// Display form of the Plotters backend error.
        message: String,
    },
}

impl fmt::Display for PlotError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Setup { source, message } => {
                write!(
                    formatter,
                    "failed to setup {source} with message: {message} ",
                )
            }

            Self::Csv { path, source } => {
                write!(formatter, "failed to read {}: {source}", path.display())
            }
            Self::InvalidData { path, message } => {
                write!(
                    formatter,
                    "invalid result data in {}: {message}",
                    path.display()
                )
            }
            Self::CreateDirectory { path, source } => write!(
                formatter,
                "failed to create plot directory {}: {source}",
                path.display()
            ),
            Self::Draw { path, message } => {
                write!(formatter, "failed to draw {}: {message}", path.display())
            }
        }
    }
}

impl Error for PlotError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Setup { .. } => None,
            Self::Csv { source, .. } => Some(source),
            Self::CreateDirectory { source, .. } => Some(source),
            Self::InvalidData { .. } | Self::Draw { .. } => None,
        }
    }
}
