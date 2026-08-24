pub mod inputs;

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
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    Parse {
        path: PathBuf,
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
pub fn default_inputs() -> SimulationInputs {
    SimulationInputs::default()
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
