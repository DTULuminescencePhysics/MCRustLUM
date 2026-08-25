// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Contains various structs to store all input data. currently this is in a single module
//! however as more functionality is added that may change. 

use common::constants::time::TimeUnit;
use common::constants::temperature::TemperatureUnit;
use common::numeric::{Float, TimeFloat};

/// Geometry, site density, and boundary settings used to construct a cube.
///
/// Lengths are expressed in metres and `density` is the volumetric trap
/// density. Hole and bandtail counts are ratios per generated trap.
#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct CubeSpecification {
    /// Unit-cell height in metres.
    pub uc_h: Float,
    /// Unit-cell width in metres.
    pub uc_w: Float,
    /// Unit-cell length in metres.
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
    /// Mean localised-state energies.
    pub e_loc: Vec<Float>,
    /// Mean conduction-band activation energies.
    pub e_cb: Vec<Float>,
    /// Standard deviation associated with each localised energy.
    pub e_loc_sigma: Vec<Float>,
    /// Standard deviation associated with each conduction-band energy.
    pub e_cb_sigma: Vec<Float>,
}

impl Default for TrapEnergies {
    fn default() -> Self {
        Self {
            e_loc: vec![1.2],
            e_cb: vec![2.0],
            e_loc_sigma: vec![0.0],
            e_cb_sigma: vec![0.0],
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
    /// Enable ground-state localised retrapping.
    pub gs_retrap: bool,
    /// Enable excited-state localised retrapping.
    pub es_retrap: bool,
    /// Enable variable-range hopping when supported by the simulation.
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
            gs_retrap: true,
            es_retrap: true,
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
    /// Enable retrapping from the conduction band.
    pub retrap: bool,
    /// Ground-state frequency factors.
    pub s_gs: Vec<Float>,
    /// Excited-state frequency factors.
    pub s_es: Vec<Float>,
    /// General-order kinetic exponents.
    pub mu: Vec<Float>,
    /// Relative retrapping strengths for each configured trap family.
    pub retrap_ratio: Vec<Float>,
}

impl Default for DeLocalisedInputs {
    fn default() -> Self {
        Self {
            gs_cb: true,
            es_cb: true,
            retrap: true,
            s_gs: vec![1.2e12],
            s_es: vec![1.2e12],
            mu: vec![0.1],
            retrap_ratio: vec![0.0],
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
    /// Allow recombination while the system is being filled.
    pub cmbn_whn_fll: bool,
    /// Recombination prefactors used during filling.
    pub recm_pre_fll: Vec<Float>,
}

impl Default for FillingInputs {
    fn default() -> Self {
        Self {
            fill: true,
            d0: vec![400.0],
            d_dot: vec![1.0],
            dd_unit: TimeUnit::Second,
            cmbn_whn_fll: false,
            recm_pre_fll: vec![0.0],
        }
    }
}


/// All input groups required to configure a simulation.
///
/// This is the top-level structure represented by an input TOML file. Missing
/// groups and missing values within a group use their corresponding defaults.
///
/// ```
/// let mut inputs = io::SimulationInputs::default();
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
    /// Localised transition configuration.
    pub localised: LocalisedInputs,
    /// Delocalised transition configuration.
    pub delocalised: DeLocalisedInputs,
    /// Dose-driven filling configuration.
    pub filling: FillingInputs,
}

impl Default for SimulationInputs {
    fn default() -> Self {
        Self {
            cube: CubeSpecification::default(),
            time_temperature: TimeTempSpecification::default(),
            trap_energies: TrapEnergies::default(),
            localised: LocalisedInputs::default(),
            delocalised: DeLocalisedInputs::default(),
            filling: FillingInputs::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            }
        );

        assert_eq!(
            LocalisedInputs::default(),
            LocalisedInputs {
                gs_tun: true,
                es_tun: true,
                gs_retrap: true,
                es_retrap: true,
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
                retrap: true,
                s_gs: vec![1.2e12],
                s_es: vec![1.2e12],
                mu: vec![0.1],
                retrap_ratio: vec![0.0],
            }
        );

        assert_eq!(
            FillingInputs::default(),
            FillingInputs {
                fill: true,
                d0: vec![400.0],
                d_dot: vec![1.0],
                dd_unit: TimeUnit::Second,
                cmbn_whn_fll: false,
                recm_pre_fll: vec![0.0],
            }
        );
    }
}
