// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Strongly typed configuration for crystal geometry and charge-transfer physics.
//!
//! Each structure corresponds to a TOML section. Vector-valued physical
//! parameters either contain one value shared by every experiment or one
//! value per experiment; the Monte Carlo setup validates that indexing when a
//! run is constructed.


use common::constants::{metric::metric_abbrev, quantities::{Dimension, Quantity, Unit}, temperature::TemperatureUnit};
use serde::{de::DeserializeOwned, Deserialize};
use std::{fs, path::Path};
type InputResult<T> = Result<T, String>;


#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SharedSettings {
    experiment_list: Vec<String>,
    cube: String,
    time_temperature: String,
    transitions: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ExperimentDefinition {
    filling: String,
    band: String,
    trap_types: Vec<String>,
    total_traps: usize,
    // Missing values inherit the shared setting.
    cube: Option<String>,
    time_temperature: Option<String>,
    transitions: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CubeDefinition {
    pub unit_cell_height: Quantity,
    pub unit_cell_width: Quantity,
    pub unit_cell_length: Quantity,
    pub cube_height: Quantity,
    pub cube_width: Quantity,
    pub cube_length: Quantity,
    pub hole_count: usize,
    pub bandtail_count: usize,
    pub periodic: bool,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeTemperatureDefinition {
    pub times: Quantity<Vec<f64>>,
    pub temperatures: Quantity<Vec<f64>>,
}

// These definitions contain only scalar quantities.
// The macro keeps their declarations short while retaining named,
// strongly typed fields.
macro_rules! quantity_entity {
    ($name:ident { $($field:ident),+ $(,)? }) => {
        #[derive(Debug, Clone, Deserialize)]
        #[serde(deny_unknown_fields)]
        pub struct $name {
            $(pub $field: Quantity,)+
        }
    };
}

quantity_entity!(TrapDefinition {
    density,
    initial_fill,
    ground_to_excited_state_energy,
    ground_to_conduction_band_energy,
    s_ground_to_excited_frequency,
    s_excited_to_ground_frequency,
    b_ground_state_to_hole_frequency,
    b_excited_state_to_hole_frequency,
    alpha_ground_state_barrier_constant,
    alpha_excited_state_barrier_constant,
    s_ground_state_to_conduction_band_frequency,
    s_excited_state_to_conduction_band_frequency,
});

quantity_entity!(BandDefinition {
    conduction_band_hole_to_trap_preference,
    valence_band_hole_to_trap_preference,
    mu_conduction_band_distance,
    mu_valence_band_distance,
});

quantity_entity!(FillingDefinition { d0, d_dot });

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TransitionDefinition {
    pub ground_state_tunnelling: bool,
    pub excited_state_tunnelling: bool,
    pub vrh: bool,
    pub ground_state_to_conduction_band: bool,
    pub excited_state_to_conduction_band: bool,
    pub filling: bool,
    pub conduction_band_retrapping: bool,
    pub ground_state_tunnelling_retrapping: bool,
    pub excited_state_tunnelling_retrapping: bool,
    pub filling_retrapping: bool,
    pub gaussian_kernel: bool,
}

#[derive(Debug, Clone)]
pub struct ResolvedExperiment {
    pub name: String,
    pub total_traps: usize,
    pub cube: CubeDefinition,
    pub time_temperature: TimeTemperatureDefinition,
    pub filling: FillingDefinition,
    pub band: BandDefinition,
    pub transitions: TransitionDefinition,

    // Preserve both the family name and its definition.
    pub traps: Vec<(String, TrapDefinition)>,
}



fn decode<T: DeserializeOwned>(
    table: &toml::Table,
    name: &str,
    context: &str,
) -> InputResult<T> {
    let value = table
        .get(name)
        .ok_or_else(|| format!("{context}: missing [{name}]"))?;

    value
        .clone()
        .try_into::<T>()
        .map_err(|error| format!("{context}: [{name}]: {error}"))
}

fn entity<T: DeserializeOwned>(
    root: &toml::Table,
    section: &str,
    name: &str,
    experiment: &str,
) -> InputResult<T> {
    let registry = root
        .get(section)
        .and_then(toml::Value::as_table)
        .ok_or_else(|| {
            format!("{experiment}: missing entity registry [{section}]")
        })?;

    decode(
        registry,
        name,
        &format!("{experiment}: reference to {section}.{name}"),
    )
}

pub fn read_experiments(
    path: impl AsRef<Path>,
    experiment_count: usize,
) -> InputResult<Vec<ResolvedExperiment>> {
    let path = path.as_ref();

    let contents = fs::read_to_string(path)
        .map_err(|error| format!("{}: {error}", path.display()))?;

    parse_experiments(&contents, experiment_count)
        .map_err(|error| format!("{}: {error}", path.display()))
}

pub fn parse_experiments(
    contents: &str,
    experiment_count: usize,
) -> InputResult<Vec<ResolvedExperiment>> {
    let root: toml::Table = toml::from_str(contents)
        .map_err(|error| format!("invalid TOML: {error}"))?;

    let shared: SharedSettings =
        decode(&root, "shared_settings", "configuration")?;

    if experiment_count == 0 {
        return Err("experiment_count must be positive".into());
    }

    if experiment_count > shared.experiment_list.len() {
        return Err(format!(
            "requested {experiment_count} experiments, \
             but experiment_list contains {}",
            shared.experiment_list.len(),
        ));
    }

    shared
        .experiment_list
        .iter()
        .take(experiment_count)
        .map(|name| {
            let definition: ExperimentDefinition =
                decode(&root, name, "experiment")?;

            if definition.total_traps == 0 {
                return Err(format!("{name}: total_traps must be positive"));
            }

            if definition.trap_types.is_empty() {
                return Err(format!("{name}: trap_types must not be empty"));
            }

            let cube_name = definition
                .cube
                .as_deref()
                .unwrap_or(&shared.cube);

            let profile_name = definition
                .time_temperature
                .as_deref()
                .unwrap_or(&shared.time_temperature);

            let transitions_name = definition
                .transitions
                .as_deref()
                .unwrap_or(&shared.transitions);

            let profile: TimeTemperatureDefinition =
                entity(&root, "time_temperature", profile_name, name)?;

            if profile.times.value.len() < 2
                || profile.times.value.len()
                    != profile.temperatures.value.len()
            {
                return Err(format!(
                    "{name}: time_temperature.{profile_name} needs \
                     at least two paired time/temperature values"
                ));
            }

            // Array distributions need separate semantics; reject them here.
            if profile.times.distribute || profile.temperatures.distribute {
                return Err(format!(
                    "{name}: distributions on profile arrays are unsupported"
                ));
            }

            let traps = definition
                .trap_types
                .iter()
                .map(|trap_name| {
                    let trap =
                        entity(&root, "traps", trap_name, name)?;
                    Ok((trap_name.clone(), trap))
                })
                .collect::<InputResult<Vec<_>>>()?;

            Ok(ResolvedExperiment {
                name: name.clone(),
                total_traps: definition.total_traps,
                cube: entity(&root, "cubes", cube_name, name)?,
                time_temperature: profile,
                filling: entity(
                    &root, "filling", &definition.filling, name,
                )?,
                band: entity(&root, "bands", &definition.band, name)?,
                transitions: entity(
                    &root, "transitions", transitions_name, name,
                )?,
                traps,
            })
        })
        .collect()
}

impl ResolvedExperiment {

    pub fn alter_time_dependent_units(&mut self, time_string: &Unit) -> Result<(), String> {
        let kind = match time_string.dimension{
            Dimension::Time { kind } => kind,
            _ => return Err(format!("Should be a time kind but instead {:?}",time_string.dimension))

        };
        self.filling.d_dot.change_unit(
            Unit {
                metric: metric_abbrev::STANDARD,
                dimension: Dimension::DoseRate{kind}}
        )?;
        
        let new_unit = Unit {
                metric: metric_abbrev::STANDARD,
                dimension: Dimension::Frequency {kind}};
        
        for (_, trap) in &mut self.traps {
            for quantity in [
                 &mut trap.s_ground_to_excited_frequency,
                &mut trap.s_excited_to_ground_frequency,
                &mut trap.b_ground_state_to_hole_frequency,
                &mut trap.b_excited_state_to_hole_frequency,
                &mut trap.s_ground_state_to_conduction_band_frequency,
                &mut trap.s_excited_state_to_conduction_band_frequency,
            ]{
                quantity.change_unit(new_unit.clone())?;
            }
        }
        
        Ok(())

    }

    pub fn alter_distance_dependent_units(&mut self, distance_string: &metric_abbrev) -> Result<(), String> {
        let new_unit = Unit { metric: distance_string.clone(), dimension: Dimension::InverseLength };
        for quantity in [
            &mut self.band.mu_conduction_band_distance,
            &mut self.band.mu_valence_band_distance
        ] {
            quantity.change_unit(new_unit.clone())?;

        }
       
        let new_unit= Unit { metric:distance_string.clone(), dimension: Dimension::InverseLength };
        for (_, trap) in &mut self.traps {
            for quantity in [
                &mut trap.density,
                &mut trap.alpha_ground_state_barrier_constant,
                &mut trap.alpha_excited_state_barrier_constant
            ] {
                quantity.change_unit(new_unit.clone())?;
            }
        }
        Ok(())
    }
}


pub fn read_inputs(
    path: impl AsRef<Path>,
    experiment_count: usize,) -> Result<Vec<ResolvedExperiment>, String> {

        let raw_exp = read_experiments(path, experiment_count)?;
        let mut experiments: Vec<ResolvedExperiment> = Vec::with_capacity(experiment_count);
        for i in 0..experiment_count {
            let mut experiment = raw_exp[i].clone();
            let time_string = &raw_exp[i].time_temperature.times.unit;
            let distance_string = &raw_exp[i].cube.cube_height.unit.metric;
            experiment.alter_time_dependent_units(time_string)?;
            experiment.alter_distance_dependent_units(distance_string)?;
            experiment.time_temperature.temperatures.change_unit(Unit{ metric: metric_abbrev::STANDARD, dimension: Dimension::Temperature{kind: TemperatureUnit::Kelvin}})?;
            experiments.push(experiment);
        }

        Ok(experiments)


}
