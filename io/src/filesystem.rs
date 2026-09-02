// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Creation of the directory layout used by a simulation run.

use std::error::Error;
use std::ffi::OsStr;
use std::fmt;
use std::fs;
use std::path::{Component, Path, PathBuf};

/// An error produced while preparing an experiment directory.
#[derive(Debug)]
pub enum FilesystemError {
    /// A supplied folder name was empty, nested, or otherwise unsafe.
    InvalidFolderName { name: PathBuf },
    /// A filesystem operation failed.
    Operation {
        action: &'static str,
        path: PathBuf,
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

fn operation_error(
    action: &'static str,
    path: &Path,
) -> impl FnOnce(std::io::Error) -> FilesystemError {
    let path = path.to_path_buf();
    move |source| FilesystemError::Operation {
        action,
        path,
        source,
    }
}

fn validate_folder_name(name: &OsStr) -> Result<(), FilesystemError> {
    let path = Path::new(name);
    let mut components = path.components();
    let valid = matches!(components.next(), Some(Component::Normal(component)) if component == name)
        && components.next().is_none();

    if valid {
        Ok(())
    } else {
        Err(FilesystemError::InvalidFolderName {
            name: path.to_path_buf(),
        })
    }
}

fn create_numbered_directory(run_directory: &Path) -> Result<PathBuf, FilesystemError> {
    for number in 1..=usize::MAX {
        let experiment_directory = run_directory.join(format!("experiment_{number}"));
        match fs::create_dir(&experiment_directory) {
            Ok(()) => return Ok(experiment_directory),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(source) => {
                return Err(FilesystemError::Operation {
                    action: "create experiment directory",
                    path: experiment_directory,
                    source,
                });
            }
        }
    }

    Err(FilesystemError::ExperimentNumberExhausted)
}

fn create_experiment_directory(
    starting_directory: &Path,
    folder_name: Option<&OsStr>,
) -> Result<PathBuf, FilesystemError> {
    let source_input = starting_directory.join("input.toml");
    let input_metadata =
        fs::metadata(&source_input).map_err(operation_error("read input file", &source_input))?;
    if !input_metadata.is_file() {
        return Err(FilesystemError::Operation {
            action: "read input file",
            path: source_input,
            source: std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "input.toml is not a regular file",
            ),
        });
    }

    let run_directory = starting_directory.join("run");
    fs::create_dir_all(&run_directory)
        .map_err(operation_error("create run directory", &run_directory))?;

    let experiment_directory = match folder_name {
        Some(name) => {
            validate_folder_name(name)?;
            let directory = run_directory.join(name);
            fs::create_dir(&directory)
                .map_err(operation_error("create experiment directory", &directory))?;
            directory
        }
        None => create_numbered_directory(&run_directory)?,
    };

    let temporary_directory = experiment_directory.join("tmp");
    if let Err(source) = fs::create_dir(&temporary_directory) {
        let _ = fs::remove_dir(&experiment_directory);
        return Err(FilesystemError::Operation {
            action: "create temporary-output directory",
            path: temporary_directory,
            source,
        });
    }

    let destination_input = experiment_directory.join("input.toml");
    if let Err(source) = fs::copy(&source_input, &destination_input) {
        let _ = fs::remove_dir(&temporary_directory);
        let _ = fs::remove_dir(&experiment_directory);
        return Err(FilesystemError::Operation {
            action: "copy input file to",
            path: destination_input,
            source,
        });
    }

    Ok(experiment_directory)
}

/// Prepare an experiment folder and make it the process working directory.
///
/// A supplied name creates `run/<name>`. With no name, the first available
/// `run/experiment_N` directory is selected, starting at one. The new
/// experiment directory contains `tmp/` and a copy of the original
/// `input.toml`.
pub fn prepare_experiment_directory(
    folder_name: Option<&OsStr>,
) -> Result<PathBuf, FilesystemError> {
    let starting_directory =
        std::env::current_dir().map_err(|source| FilesystemError::Operation {
            action: "determine current directory",
            path: PathBuf::from("."),
            source,
        })?;
    let experiment_directory = create_experiment_directory(&starting_directory, folder_name)?;

    std::env::set_current_dir(&experiment_directory).map_err(operation_error(
        "change working directory to",
        &experiment_directory,
    ))?;

    Ok(experiment_directory)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_starting_directory(label: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should follow the Unix epoch")
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "mcrustlum_filesystem_{label}_{}_{}",
            std::process::id(),
            unique,
        ));
        fs::create_dir(&directory).unwrap();
        fs::write(directory.join("input.toml"), "[cube]\nx = 1.0\n").unwrap();
        directory
    }

    #[test]
    fn named_experiment_contains_tmp_and_an_input_copy() {
        let starting_directory = temporary_starting_directory("named");

        let experiment_directory =
            create_experiment_directory(&starting_directory, Some(OsStr::new("my_experiment")))
                .unwrap();

        assert_eq!(
            experiment_directory,
            starting_directory.join("run/my_experiment")
        );
        assert!(experiment_directory.join("tmp").is_dir());
        assert_eq!(
            fs::read_to_string(experiment_directory.join("input.toml")).unwrap(),
            "[cube]\nx = 1.0\n"
        );
        fs::remove_dir_all(starting_directory).unwrap();
    }

    #[test]
    fn unnamed_experiments_use_the_first_available_number() {
        let starting_directory = temporary_starting_directory("numbered");
        fs::create_dir(starting_directory.join("run")).unwrap();
        fs::create_dir(starting_directory.join("run/experiment_1")).unwrap();

        let experiment_directory = create_experiment_directory(&starting_directory, None).unwrap();

        assert_eq!(
            experiment_directory,
            starting_directory.join("run/experiment_2")
        );
        fs::remove_dir_all(starting_directory).unwrap();
    }

    #[test]
    fn rejects_nested_names_and_existing_named_experiments() {
        let starting_directory = temporary_starting_directory("invalid");

        let invalid =
            create_experiment_directory(&starting_directory, Some(OsStr::new("parent/child")))
                .unwrap_err();
        assert!(matches!(invalid, FilesystemError::InvalidFolderName { .. }));

        create_experiment_directory(&starting_directory, Some(OsStr::new("existing"))).unwrap();
        let existing =
            create_experiment_directory(&starting_directory, Some(OsStr::new("existing")))
                .unwrap_err();
        assert!(matches!(existing, FilesystemError::Operation { .. }));

        fs::remove_dir_all(starting_directory).unwrap();
    }

    #[test]
    fn missing_input_does_not_consume_an_experiment_number() {
        let starting_directory = temporary_starting_directory("missing_input");
        fs::remove_file(starting_directory.join("input.toml")).unwrap();

        let error = create_experiment_directory(&starting_directory, None).unwrap_err();

        assert!(matches!(error, FilesystemError::Operation { .. }));
        assert!(!starting_directory.join("run").exists());
        fs::remove_dir(starting_directory).unwrap();
    }
}
