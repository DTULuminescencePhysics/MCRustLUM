use common::constants::time::{TemperatureUnit, TimeUnit};
use common::numeric::{Float, TimeFloat};

#[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct CubeSpecification {
    pub uc_h: Float,
    pub uc_w: Float,
    pub uc_l: Float,
    pub x: Float,
    pub y: Float,
    pub z: Float,
    pub density: Float,
    pub hole_count: usize,
    pub bandtail_count: usize,
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

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct TimeTempSpecification {
    pub times: Vec<TimeFloat>,
    pub temperatures: Vec<Float>,
    pub time_unit: TimeUnit,
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

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct TrapEnergies {
    pub e_loc: Vec<Float>,
    pub e_cb: Vec<Float>,
    pub e_loc_sigma: Vec<Float>,
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

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct LocalisedInputs {
    pub gs_tun: bool,
    pub es_tun: bool,
    pub retrap: bool,
    pub vrh: bool,
    pub b_gs: Vec<Float>,
    pub b_es: Vec<Float>,
    pub alpha_gs: Vec<Float>,
    pub alpha_es: Vec<Float>,
}

impl Default for LocalisedInputs {
    fn default() -> Self {
        Self {
            gs_tun: true,
            es_tun: true,
            retrap: true,
            vrh: false,
            b_gs: vec![1.2e12],
            b_es: vec![1.2e12],
            alpha_gs: vec![9.0e12],
            alpha_es: vec![9.0e9],
        }
    }
}

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct DeLocalisedInputs {
    pub gs_cb: bool,
    pub es_cb: bool,
    pub s_gs: Vec<Float>,
    pub s_es: Vec<Float>,
    pub mu: Vec<Float>,
    pub retrap_ratio: Vec<Float>,
}

impl Default for DeLocalisedInputs {
    fn default() -> Self {
        Self {
            gs_cb: true,
            es_cb: true,
            s_gs: vec![1.2e12],
            s_es: vec![1.2e12],
            mu: vec![0.1],
            retrap_ratio: vec![0.0],
        }
    }
}

#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct FillingInputs {
    pub fill: bool,
    pub d0: Vec<Float>,
    pub d_dot: Vec<Float>,
    pub dd_unit: TimeUnit,
    pub cmbn_whn_fll: bool,
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
#[derive(Debug, Clone, PartialEq, serde::Deserialize)]
#[serde(default)]
pub struct SimulationInputs {
    pub cube: CubeSpecification,
    pub time_temperature: TimeTempSpecification,
    pub trap_energies: TrapEnergies,
    pub localised: LocalisedInputs,
    pub delocalised: DeLocalisedInputs,
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
                retrap: true,
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
