// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Strongly typed configuration for crystal geometry and charge-transfer physics.
//!
//! Named TOML entities describe geometry, trap families, capture, filling, and
//! temperature profiles. [`crate::inputs::read_experiments`] resolves references and shared
//! defaults; [`crate::inputs::read_inputs`] also normalizes the units needed by simulation.
//! A [`crate::inputs::ResolvedExperiment`] contains concrete definitions in configured order.


use common::constants::{metric::MetricAbbrev, quantities::{Dimension, Quantity, Unit}, temperature::TemperatureUnit};
use serde::{de::DeserializeOwned, Deserialize};
use std::{fs, path::Path};
type InputResult<T> = Result<T, String>;


#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
/// Ordered experiment names and default references to shared geometry, profiles, and transitions.
struct SharedSettings {
    experiment_list: Vec<String>,
    cube: String,
    time_temperature: String,
    transitions: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
/// Unresolved experiment references read from a named TOML section.
///
/// Optional cube, profile, and transition names override shared defaults.
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
/// Unit-bearing cell/crystal dimensions and site counts for one geometry.
///
/// `periodic` selects periodic boundaries. In density-based construction, hole
/// and band-tail counts are multipliers per trap; the resolved experiment
/// supplies the minimum total trap population.
pub struct CubeDefinition {
    /// Height of one unit cell, with an explicit length unit.
    pub unit_cell_height: Quantity,
    /// Width of one unit cell, with an explicit length unit.
    pub unit_cell_width: Quantity,
    /// Length of one unit cell, with an explicit length unit.
    pub unit_cell_length: Quantity,
    /// Crystal height; its metric prefix sets runtime spatial parameter units.
    pub cube_height: Quantity,
    /// Crystal width, with an explicit length unit.
    pub cube_width: Quantity,
    /// Crystal length, with an explicit length unit.
    pub cube_length: Quantity,
    /// Hole sites per trap in the density-based geometry constructor.
    pub hole_count: usize,
    /// Band-tail sites per trap in the density-based geometry constructor.
    pub bandtail_count: usize,
    /// Whether distances wrap across opposite crystal faces.
    pub periodic: bool,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
/// Paired time and temperature quantities describing a piecewise-linear profile.
///
/// Arrays carry their own units, must have equal lengths, and need at least two
/// points. Distributions on these arrays are rejected by the experiment parser.
pub struct TimeTemperatureDefinition {
    /// Time coordinates in the profile time unit, paired with temperature coordinates.
    pub times: Quantity<Vec<f64>>,
    /// Temperature coordinates paired with times, in the accompanying temperature unit.
    pub temperatures: Quantity<Vec<f64>>,
}

// These definitions contain only scalar quantities.
// The macro keeps their declarations short while retaining named,
// strongly typed fields.
macro_rules! quantity_entity {
    ($(#[$type_doc:meta])* $name:ident { $($(#[$field_doc:meta])* $field:ident),+ $(,)? }) => {
        #[derive(Debug, Clone, Deserialize)]
        #[serde(deny_unknown_fields)]
        $(#[$type_doc])*
        pub struct $name {
            $(
                $(#[$field_doc])*
                pub $field: Quantity,
            )+
        }
    };
}

quantity_entity!(
    /// Physical parameters for one named trap family, each retaining its unit and distribution.
    TrapDefinition {
        /// Relative/physical trap-family density used when apportioning sites; retains its input unit.
        density,
        /// Initial occupied fraction for this family (dimensionless).
        initial_fill,
        /// Energy separation between localized ground and excited states.
        ground_to_excited_state_energy,
        /// Ground-state activation energy for release to the conduction band.
        ground_to_conduction_band_energy,
        /// Attempt frequency for local ground-to-excited thermal transitions.
        s_ground_to_excited_frequency,
        /// Relaxation frequency from the excited state to the ground state.
        s_excited_to_ground_frequency,
        /// Ground-state tunnelling prefactor for transfer to hole sites.
        b_ground_state_to_hole_frequency,
        /// Excited-state tunnelling prefactor for transfer to hole sites.
        b_excited_state_to_hole_frequency,
        /// Ground-state tunnelling decay constant in inverse length units.
        alpha_ground_state_barrier_constant,
        /// Excited-state tunnelling decay constant in inverse length units.
        alpha_excited_state_barrier_constant,
        /// Ground-state attempt frequency for thermal release to the conduction band.
        s_ground_state_to_conduction_band_frequency,
        /// Excited-state attempt frequency for thermal release to the conduction band.
        s_excited_state_to_conduction_band_frequency,
});

quantity_entity!(
    /// Band-capture preference ratios and characteristic distances for one named configuration.
    BandDefinition {
        /// Relative hole/trap capture preference for conduction-band electrons.
        conduction_band_hole_to_trap_preference,
        /// Relative hole/trap capture preference for the filling capture model.
        valence_band_hole_to_trap_preference,
        /// Characteristic length for conduction-band capture distance weighting.
        mu_conduction_band_distance,
        /// Characteristic length for filling capture distance weighting.
        mu_valence_band_distance,
});

quantity_entity!(
    /// Characteristic saturation dose and applied dose rate for irradiation filling.
    FillingDefinition {
        /// Characteristic dose scale in the filling rate equation.
        d0,
        /// Applied dose per unit time.
        d_dot,
});

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
/// Flags enabling microscopic charge-transfer pathways and capture kernels.
///
/// Ground/excited tunnelling, thermal release, filling, and their retrapping
/// branches are selected independently. `gaussian_kernel` selects distance weighting.
pub struct TransitionDefinition {
    /// Enable localised tunnelling from trap ground states.
    pub ground_state_tunnelling: bool,
    /// Enable localised tunnelling from trap excited states.
    pub excited_state_tunnelling: bool,
    /// Variable-range hopping selection flag supplied to transition construction.
    pub vrh: bool,
    /// Enable thermal release from ground states to the conduction band.
    pub ground_state_to_conduction_band: bool,
    /// Enable thermal release from excited states to the conduction band.
    pub excited_state_to_conduction_band: bool,
    /// Irradiation/filling configuration, or its enabling flag in a transition selection.
    pub filling: bool,
    /// Allow conduction-band electrons to be recaptured by empty traps.
    pub conduction_band_retrapping: bool,
    /// Allow localised ground-state transfer into an empty trap.
    pub ground_state_tunnelling_retrapping: bool,
    /// Allow localised excited-state transfer into an empty trap.
    pub excited_state_tunnelling_retrapping: bool,
    /// Enable the retrapping branch of filling-time capture.
    pub filling_retrapping: bool,
    /// Select Gaussian distance weighting for capture competition.
    pub gaussian_kernel: bool,
}

#[derive(Debug, Clone)]
/// One experiment with every named configuration reference resolved.
///
/// Trap family order is preserved. `read_experiments` retains supplied units;
/// `read_inputs` additionally normalizes time-dependent parameters, spatial
/// parameters, and profile temperatures for simulation.
pub struct ResolvedExperiment {
    /// Experiment name from the ordered shared experiment list.
    pub name: String,
    /// Total trap population apportioned among the configured families.
    pub total_traps: usize,
    /// Resolved geometry and site-count configuration.
    pub cube: CubeDefinition,
    /// Resolved time-temperature control points and their units.
    pub time_temperature: TimeTemperatureDefinition,
    /// Irradiation/filling configuration, or its enabling flag in a transition selection.
    pub filling: FillingDefinition,
    /// Capture preference ratios and characteristic distances.
    pub band: BandDefinition,
    /// Flags selecting the enabled physical pathways.
    pub transitions: TransitionDefinition,

    // Preserve both the family name and its definition.
    /// Ordered trap family names paired with their physical definitions.
    pub traps: Vec<(String, TrapDefinition)>,
}



/// Deserialize a named table entry, adding section and caller context to errors.
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

/// Resolve a named entity within a registry and deserialize it for an experiment.
///
/// Missing registries, missing names, and invalid entity values return contextual errors.
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

/// Read a TOML file and resolve the first `experiment_count` configured experiments.
///
/// Preserves input units. File, syntax, reference, and validation errors include
/// the source path. Use `read_inputs` for simulation-ready unit normalization.
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

/// Resolve experiments from TOML text in `shared_settings.experiment_list` order.
///
/// Uses per-experiment overrides for cube, profile, and transitions, otherwise
/// inheriting shared references. Rejects zero/excessive experiment counts,
/// empty trap families, zero trap populations, unpaired/short profiles, and
/// profile-array distributions. Unit normalization is performed separately.
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

    /// Convert dose rate and trap frequencies to the profile time unit in place.
    ///
    /// Requires a time dimension; invalid dimensions or incompatible conversions
    /// return an error. Earlier conversions remain applied if a later one fails.
    pub fn alter_time_dependent_units(&mut self, time_string: &Unit) -> Result<(), String> {
        let kind = match time_string.dimension{
            Dimension::Time { kind } => kind,
            _ => return Err(format!("Should be a time kind but instead {:?}",time_string.dimension))

        };
        self.filling.d_dot.change_unit(
            Unit {
                metric: MetricAbbrev::STANDARD,
                dimension: Dimension::DoseRate{kind}}
        )?;
        
        let new_unit = Unit {
                metric: MetricAbbrev::STANDARD,
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

    /// Convert capture distances, trap densities, and decay constants to a common prefix.
    ///
    /// Uses `distance_string` for length, inverse volume, and inverse length units.
    /// Conversion errors propagate; earlier updates are retained on failure.
    pub fn alter_distance_dependent_units(&mut self, distance_string: &MetricAbbrev) -> Result<(), String> {
        let new_unit = Unit { metric: distance_string.clone(), dimension: Dimension::Length };
        for quantity in [
            &mut self.band.mu_conduction_band_distance,
            &mut self.band.mu_valence_band_distance
        ] {
            quantity.change_unit(new_unit.clone())?;

        }
        let new_unit= Unit { metric: distance_string.clone(), dimension: Dimension::Density };
        
        for (_, trap) in &mut self.traps {
            for quantity in [
                &mut trap.density,
            ] {
                quantity.change_unit(new_unit.clone())?;
            }
        }
        let new_unit= Unit { metric:distance_string.clone(), dimension: Dimension::InverseLength };
        for (_, trap) in &mut self.traps {
            for quantity in [
                &mut trap.alpha_ground_state_barrier_constant,
                &mut trap.alpha_excited_state_barrier_constant
            ] {
                quantity.change_unit(new_unit.clone())?;
            }
        }
        Ok(())
    }
}


/// Load resolved experiments and normalize their physical units for simulation.
///
/// Converts rates to each profile's time unit, spatial parameters to the cube
/// height's metric prefix, and profile temperatures to Kelvin. Returns errors
/// from reading, resolving, or converting the requested configurations.
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
            experiment.time_temperature.temperatures.change_unit(Unit{ metric: MetricAbbrev::STANDARD, dimension: Dimension::Temperature{kind: TemperatureUnit::Kelvin}})?;
            experiments.push(experiment);
        }

        Ok(experiments)


}
