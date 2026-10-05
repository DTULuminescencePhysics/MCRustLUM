// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Strongly typed configuration for crystal geometry and charge-transfer physics.
//!
//! Each structure corresponds to a TOML section. Vector-valued physical
//! parameters either contain one value shared by every experiment or one
//! value per experiment; the Monte Carlo setup validates that indexing when a
//! run is constructed.

use crate::errors::InputError;
use common::constants::temperature::TemperatureUnit;
use common::constants::time::TimeUnit;
use common::numeric::{Float, TimeFloat};
use std::fs::read_to_string;
use std::path::Path;

/// Geometry, site density, and boundary settings used to construct a cube.
///
/// Lengths are expressed in metres and `density` is the volumetric trap
/// density. Hole and bandtail counts are ratios per generated trap.
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct CubeSpecification {
    /// Unit-cell height in metres, retained for future lattice-aware placement.
    ///
    /// Current site generation samples a continuum and does not use this value.
    pub uc_h: Float,
    /// Unit-cell width in metres, retained for future lattice-aware placement.
    ///
    /// Current site generation samples a continuum and does not use this value.
    pub uc_w: Float,
    /// Unit-cell length in metres, retained for future lattice-aware placement.
    ///
    /// Current site generation samples a continuum and does not use this value.
    pub uc_l: Float,
    /// Cube extent along the x axis in metres.
    pub x: Float,
    /// Cube extent along the y axis in metres.
    pub y: Float,
    /// Cube extent along the z axis in metres.
    pub z: Float,
    /// Number of electron traps per cubic metre.
    pub density: Float,
    /// Number of holes generated per electron trap.
    pub hole_count: usize,
    /// Number of bandtail states generated per electron trap.
    pub bandtail_count: usize,
    /// Whether distances wrap across opposite cube faces.
    pub periodic: bool,
}

impl Default for CubeSpecification {
    fn default() -> Self {
        Self {
            uc_h: 1.0e-10,
            uc_w: 1.0e-10,
            uc_l: 1.0e-10,
            x: 7.5e-9,
            y: 7.5e-9,
            z: 7.5e-9,
            density: 5.22e25,
            hole_count: 1,
            bandtail_count: 0,
            periodic: true,
        }
    }
}

/// Control points and units for a piecewise-linear temperature history.
///
/// `times` and `temperatures` are paired by index and must have equal lengths.
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct TimeTempSpecification {
    /// Time coordinate for every profile control point.
    pub times: Vec<TimeFloat>,
    /// Temperature at every profile control point.
    pub temperatures: Vec<Float>,
    /// Unit applied to every value in [`Self::times`].
    pub time_unit: TimeUnit,
    /// Unit applied to every value in [`Self::temperatures`].
    pub temp_unit: TemperatureUnit,
}

impl Default for TimeTempSpecification {
    fn default() -> Self {
        Self {
            times: vec![0.0, 160.0],
            temperatures: vec![0.0, 800.0],
            time_unit: TimeUnit::Second,
            temp_unit: TemperatureUnit::Celsius,
        }
    }
}

/// Energy distributions for localised and conduction-band transitions.
///
/// Energies and their standard deviations are expressed in electronvolts.
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct TrapEnergies {
    /// Ground-to-excited localised-state energy gaps, in eV.
    pub e_loc: Vec<Float>,
    /// Ground-state activation energies for conduction-band release, in eV.
    pub e_cb: Vec<Float>,
    /// Standard deviation associated with each localised energy, in eV.
    ///
    /// This is parsed and stored but the current homogeneous layout does not
    /// yet sample an energy distribution from it.
    pub e_loc_sigma: Vec<Float>,
    /// Standard deviation associated with each conduction-band energy, in eV.
    ///
    /// This is parsed and stored but the current homogeneous layout does not
    /// yet sample an energy distribution from it.
    pub e_cb_sigma: Vec<Float>,
    /// Attempt frequency for thermal excitation from ground to excited state, in s⁻¹.
    pub s_frequency_e: Vec<Float>,
    /// Relaxation frequency from excited to ground state, in s⁻¹.
    pub s_frequency_g: Vec<Float>,
}

impl Default for TrapEnergies {
    fn default() -> Self {
        Self {
            e_loc: vec![1.2],
            e_cb: vec![2.0],
            e_loc_sigma: vec![0.0],
            e_cb_sigma: vec![0.0],
            s_frequency_e: vec![9.0e10],
            s_frequency_g: vec![9.0e10],
        }
    }
}
/// Initial trap/hole availability.
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct InitialConditions {
    /// Initial occupied-trap fraction, expressed from `0.0` through `1.0`.
    pub trap_available: Vec<Float>,
    /// Initial active-hole fraction, expressed from `0.0` through `1.0`.
    pub hole_available: Vec<Float>,
}

impl Default for InitialConditions {
    fn default() -> Self {
        Self {
            trap_available: vec![1.0],
            hole_available: vec![1.0],
        }
    }
}
/// Selection and parameters for localised tunnelling transitions.
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct LocalisedInputs {
    /// Enable ground-state tunnelling recombination.
    pub gs_tun: bool,
    /// Enable excited-state tunnelling recombination.
    pub es_tun: bool,
    /// Request variable-range hopping in a future transport model.
    ///
    /// The current standard runner stores but does not act on this flag.
    pub vrh: bool,
    /// Ground-state attempt frequencies.
    pub b_gs: Vec<Float>,
    /// Excited-state attempt frequencies.
    pub b_es: Vec<Float>,
    /// Ground-state spatial decay constants.
    pub alpha_gs: Vec<Float>,
    /// Excited-state spatial decay constants.
    pub alpha_es: Vec<Float>,
}

impl Default for LocalisedInputs {
    fn default() -> Self {
        Self {
            gs_tun: true,
            es_tun: true,
            vrh: false,
            b_gs: vec![1.2e12],
            b_es: vec![1.2e12],
            alpha_gs: vec![9.0e12],
            alpha_es: vec![9.0e9],
        }
    }
}

/// Selection and parameters for conduction-band transitions.
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct DeLocalisedInputs {
    /// Enable ground-state conduction-band release.
    pub gs_cb: bool,
    /// Enable excited-state conduction-band release.
    pub es_cb: bool,
    /// Ground-state frequency factors.
    pub s_gs: Vec<Float>,
    /// Excited-state frequency factors.
    pub s_es: Vec<Float>,
}

impl Default for DeLocalisedInputs {
    fn default() -> Self {
        Self {
            gs_cb: true,
            es_cb: true,
            s_gs: vec![1.2e12],
            s_es: vec![1.2e12],
        }
    }
}

/// Dose-driven trap-filling configuration.
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct FillingInputs {
    /// Enable dose-driven filling.
    pub fill: bool,
    /// Characteristic doses for the configured trap families.
    pub d0: Vec<Float>,
    /// Applied dose rates.
    pub d_dot: Vec<Float>,
    /// Time denominator used by [`Self::d_dot`].
    pub dd_unit: TimeUnit,
}

impl Default for FillingInputs {
    fn default() -> Self {
        Self {
            fill: true,
            d0: vec![400.0],
            d_dot: vec![1.0],
            dd_unit: TimeUnit::Second,
        }
    }
}
/// Re-trapping configuration
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct ReTrapping {
    /// Enable retrapping for electrons in the conduction band
    pub delocalised: bool,
    /// Enable ground-state localised retrapping.
    pub localised_gs: bool,
    /// Enable excited-state localised retrapping.
    pub localised_es: bool,
    /// Enable retrapping for electrons freed by the filling process
    pub filling: bool,
    /// Use the gaussian kernel for retrapping
    pub gaussian: bool,
    /// Once an electron is promoted to the conduction band what is the
    /// preference of choosing a hole over a trap
    /// 0.0 means a trap is always chosen, 1.0 means an equal likelihood
    pub cb_hole_to_trap: Vec<Float>,
    /// When a hole is produced in the valence band what is the
    /// preference of choosing to annihilate a trapped electron over localising the hole
    /// 0.0 means a hole is always chosen, 1.0 means an equal likelihood
    pub vb_trap_to_hole: Vec<Float>,
    /// Characteristic length in the distance-dependent conduction-band capture model.
    /// Values use the same length unit as the generated site coordinates.
    /// If gaussian retrapping is turned on mu_cb controls the distance the electron can travel in
    /// the conduction band  
    pub cb_mu: Vec<Float>,
    /// If gaussian retrapping is turned on mu_vb controls the distance the hole can travel in
    /// the valence band  
    pub vb_mu: Vec<Float>,
}

impl Default for ReTrapping {
    fn default() -> Self {
        Self {
            delocalised: true,
            localised_gs: true,
            localised_es: true,
            filling: true,
            gaussian: true,
            cb_hole_to_trap: vec![1.0],
            vb_trap_to_hole: vec![1.0],
            cb_mu: vec![0.1],
            vb_mu: vec![0.1],
        }
    }
}

/// All input groups required to configure a simulation.
///
/// This is the top-level structure represented by an input TOML file. Missing
/// groups and missing values within a group use their corresponding defaults.
///
/// ```
/// let mut inputs = io::inputs::SimulationInputs::default();
/// inputs.time_temperature.times = vec![0.0, 60.0];
/// inputs.time_temperature.temperatures = vec![293.15, 373.15];
/// ```
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct SimulationInputs {
    /// Crystal geometry and site-generation settings.
    pub cube: CubeSpecification,
    /// Time and temperature profile.
    pub time_temperature: TimeTempSpecification,
    /// Trap energy distributions.
    pub trap_energies: TrapEnergies,
    /// Initial trap/hole availability.
    pub initial_conditions: InitialConditions,
    /// Localised transition configuration.
    pub localised: LocalisedInputs,
    /// Delocalised transition configuration.
    pub delocalised: DeLocalisedInputs,
    /// Dose-driven filling configuration.
    pub filling: FillingInputs,
    /// ReTrapping configuration
    pub retrapping: ReTrapping,
}

impl Default for SimulationInputs {
    fn default() -> Self {
        Self {
            cube: CubeSpecification::default(),
            time_temperature: TimeTempSpecification::default(),
            trap_energies: TrapEnergies::default(),
            initial_conditions: InitialConditions::default(),
            localised: LocalisedInputs::default(),
            delocalised: DeLocalisedInputs::default(),
            filling: FillingInputs::default(),
            retrapping: ReTrapping::default(),
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
/// # fn main() -> Result<(), io::errors::InputError> {
/// let inputs = io::inputs::read_inputs("input.toml")?;
/// println!("trap density: {}", inputs.cube.density);
/// # Ok(())
/// # }
/// ```
pub fn read_inputs(path: impl AsRef<Path>) -> Result<SimulationInputs, InputError> {
    let path = path.as_ref();
    let contents = read_to_string(path).map_err(|source| InputError::Read {
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
/// let mut inputs = io::inputs::default_inputs();
/// inputs.cube.periodic = false;
/// assert!(!inputs.cube.periodic);
/// ```
pub fn default_inputs() -> SimulationInputs {
    SimulationInputs::default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn input_defaults_match_standard_configuration() {
        assert_eq!(
            CubeSpecification::default(),
            CubeSpecification {
                uc_h: 1.0e-10,
                uc_w: 1.0e-10,
                uc_l: 1.0e-10,
                x: 7.5e-9,
                y: 7.5e-9,
                z: 7.5e-9,
                density: 5.22e25,
                hole_count: 1,
                bandtail_count: 0,
                periodic: true,
            }
        );

        assert_eq!(
            TimeTempSpecification::default(),
            TimeTempSpecification {
                times: vec![0.0, 160.0],
                temperatures: vec![0.0, 800.0],
                time_unit: TimeUnit::Second,
                temp_unit: TemperatureUnit::Celsius,
            }
        );

        assert_eq!(
            TrapEnergies::default(),
            TrapEnergies {
                e_loc: vec![1.2],
                e_cb: vec![2.0],
                e_loc_sigma: vec![0.0],
                e_cb_sigma: vec![0.0],
                s_frequency_e: vec![9.0e10],
                s_frequency_g: vec![9.0e10],
            }
        );

        assert_eq!(
            LocalisedInputs::default(),
            LocalisedInputs {
                gs_tun: true,
                es_tun: true,
                vrh: false,
                b_gs: vec![1.2e12],
                b_es: vec![1.2e12],
                alpha_gs: vec![9.0e12],
                alpha_es: vec![9.0e9],
            }
        );

        assert_eq!(
            DeLocalisedInputs::default(),
            DeLocalisedInputs {
                gs_cb: true,
                es_cb: true,
                s_gs: vec![1.2e12],
                s_es: vec![1.2e12],
            }
        );

        assert_eq!(
            FillingInputs::default(),
            FillingInputs {
                fill: true,
                d0: vec![400.0],
                d_dot: vec![1.0],
                dd_unit: TimeUnit::Second,
            }
        );

        assert_eq!(
            InitialConditions::default(),
            InitialConditions {
                trap_available: vec![1.0],
                hole_available: vec![1.0],
            }
        );
        assert_eq!(
            ReTrapping::default(),
            ReTrapping {
                delocalised: true,
                localised_gs: true,
                localised_es: true,
                filling: true,
                gaussian: true,
                cb_hole_to_trap: vec![1.0],
                vb_trap_to_hole: vec![1.0],
                cb_mu: vec![0.1],
                vb_mu: vec![0.1],
            }
        );
    }

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
            InputError::Read {
                path: error_path,
                source,
            } => {
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
            InputError::Parse {
                path: error_path, ..
            } => assert_eq!(error_path, &path),
            InputError::Read { .. } => panic!("malformed TOML should produce a parse error"),
        }
        assert!(error.to_string().contains("failed to parse"));
        assert!(error.source().is_some());
    }
}
