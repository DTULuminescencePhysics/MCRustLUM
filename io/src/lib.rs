// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Loading and representing simulation configuration.
//!
//! Use [`read_inputs`] when a TOML file is supplied and [`default_inputs`]
//! when the built-in configuration is sufficient. Both paths return the same
//! [`SimulationInputs`] type, so downstream simulation code does not need to
//! know where the values came from.

/// Typed groups corresponding to the sections of an input TOML file.
pub mod inputs;

pub mod outputs;

use std::error::Error;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

pub use inputs::{
    CubeSpecification, DeLocalisedInputs, FillingInputs, LocalisedInputs, SimulationInputs,
    TimeTempSpecification, TrapEnergies,
};

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
    /// File contents were not valid TOML for [`SimulationInputs`].
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

/// Read a TOML input file and construct all grouped simulation inputs.
///
/// Values omitted from the file are filled from the `Default` implementation
/// of the relevant input structure.
///
/// # Example
///
/// ```no_run
/// # fn main() -> Result<(), io::InputError> {
/// let inputs = io::read_inputs("input.toml")?;
/// println!("trap density: {}", inputs.cube.density);
/// # Ok(())
/// # }
/// ```
pub fn read_inputs(path: impl AsRef<Path>) -> Result<SimulationInputs, InputError> {
    let path = path.as_ref();
    let contents = fs::read_to_string(path).map_err(|source| InputError::Read {
        path: path.to_path_buf(),
        source,
    })?;

    toml::from_str(&contents).map_err(|source| InputError::Parse {
        path: path.to_path_buf(),
        source,
    })
}

/// Construct a complete input set without reading a file.
///
/// ```
/// let mut inputs = io::default_inputs();
/// inputs.cube.periodic = false;
/// assert!(!inputs.cube.periodic);
/// ```
pub fn default_inputs() -> SimulationInputs {
    SimulationInputs::default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_input_path(label: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "mcrustlum_{label}_{}_{}.toml",
            std::process::id(),
            unique,
        ))
    }

    #[test]
    fn reads_the_workspace_input_file() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../input.toml");
        read_inputs(path).expect("the workspace input file should be valid");
    }

    #[test]
    fn omitted_values_use_struct_defaults() {
        let inputs: SimulationInputs = toml::from_str(
            r#"
                [cube]
                x = 1.5e-8
            "#,
        )
        .expect("a partial input file should be valid");

        assert_eq!(inputs.cube.x, 1.5e-8);
        assert_eq!(inputs.cube.y, CubeSpecification::default().y);
        assert_eq!(inputs.time_temperature, TimeTempSpecification::default());
        assert_eq!(inputs.trap_energies, TrapEnergies::default());
        assert_eq!(inputs.localised, LocalisedInputs::default());
        assert_eq!(inputs.delocalised, DeLocalisedInputs::default());
        assert_eq!(inputs.filling, FillingInputs::default());
    }

    #[test]
    fn missing_input_file_reports_the_path_and_io_source() {
        let path = temporary_input_path("missing");
        let error = read_inputs(&path).unwrap_err();

        match &error {
            InputError::Read { path: error_path, source } => {
                assert_eq!(error_path, &path);
                assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
            }
            InputError::Parse { .. } => panic!("missing file should produce a read error"),
        }
        assert!(error.to_string().contains(path.to_string_lossy().as_ref()));
        assert!(error.source().is_some());
    }

    #[test]
    fn malformed_input_file_reports_the_path_and_parse_source() {
        let path = temporary_input_path("malformed");
        fs::write(&path, "[cube\nx = not-a-number")
            .expect("temporary malformed input should be writable");

        let error = read_inputs(&path).unwrap_err();
        fs::remove_file(&path).expect("temporary input should be removable");

        match &error {
            InputError::Parse { path: error_path, .. } => assert_eq!(error_path, &path),
            InputError::Read { .. } => panic!("malformed TOML should produce a parse error"),
        }
        assert!(error.to_string().contains("failed to parse"));
        assert!(error.source().is_some());
    }
}
