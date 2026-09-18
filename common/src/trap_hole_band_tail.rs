// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! [`crate::trap_hole_band_tail::Coord`] contains the spatial coordinates of
//! each site that can contain an electron. These are then stored in vectors
//! in [`crate::trap_hole_band_tail::ElectronPlaces`]


//! Constructors validate geometry before allocating or randomly generating
//! site coordinates.
use crate::numeric::{Float, Numeric};
use crate::place_ids::{PlaceId, PlaceAvailability};
use crate::rate_equation_inputs::{
    DelocalisedTransitionInputs, FillingTransitionInputs, LocalisedTransitionInputs,
};
use rand::Rng;
/// A three-dimensional site position.
///
/// Coordinates can be supplied directly or generated within positive x, y,
/// and z limits.
#[derive(Debug, Clone, Copy)]
pub struct Coord {
    /// Position along the x axis.
    pub x: Float,
    /// Position along the y axis.
    pub y: Float,
    /// Position along the z axis.
    pub z: Float,
}

impl Coord {
    /// Create a finite, non-negative coordinate.
    pub fn new<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z) -> Result<Self, String> {
        let x = x.to_float();
        let y = y.to_float();
        let z = z.to_float();

        let coordinate = Self { x, y, z };
        coordinate.validate()?;
        Ok(coordinate)
    }

    /// Generate a coordinate in the inclusive ranges `[0, x]`, `[0, y]`, and `[0, z]`.
    pub fn random_in<X: Numeric, Y: Numeric, Z: Numeric,>(
        x: X,
        y: Y,
        z: Z,
        rng: &mut impl Rng,
    ) -> Result<Self, String> {
        let x = x.to_float();
        let y = y.to_float();
        let z = z.to_float();

        Coord::validate_dimension("x", x)?;
        Coord::validate_dimension("y", y)?;
        Coord::validate_dimension("z", z)?;

        Ok(Coord {
            x: Float::random_in(x, &mut *rng),
            y: Float::random_in(y, &mut *rng),
            z: Float::random_in(z, &mut *rng),
        })
    }

    /// Validate that one boundary dimension is finite and greater than zero.
    fn validate_dimension(name: &str, value: Float) -> Result<(), String> {
        if !value.is_finite() {
            return Err(format!("boundary {name} must be finite, got {value}"));
        }
        if value <= 0.0 {
            return Err(format!("boundary {name} must be greater than zero, got {value}"));
        }
        Ok(())
    }



    /// Validate that the coordinate is finite and non-negative.
    pub fn validate(&self) -> Result<(), String> {
        for (name, value) in [("x", self.x), ("y", self.y), ("z", self.z)] {
            if !value.is_finite() {
                return Err(format!("coordinate {name} must be finite, got {value}"));
            }
            if value < 0.0 {
                return Err(format!("coordinate {name} cannot be negative, got {value}"));
            }
        }
        Ok(())
    }

    /// Calculate the ordinary Euclidean distance to another coordinate.
    pub fn distance(&self, other: &Coord) -> Float {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2))
            .sqrt()
    }
}

/// Coordinates of the traps, holes, and bandtail states in a crystal.
#[derive(Debug, Clone)]
pub enum ElectronPlaces {
    /// Contains just traps and holes
    Standard {
        /// Positions of localised electron traps.
        traps: Vec<Coord>,
        /// Positions of hole recombination sites.
        holes: Vec<Coord>,
    },
    /// Contains traps, holes and bandtail states
    WithBandtail {
        /// Positions of localised electron traps.
        traps: Vec<Coord>,
        /// Positions of hole recombination sites.
        holes: Vec<Coord>,
        /// Positions of shallow bandtail states.
        bandtails: Vec<Coord>,
    },
}

impl ElectronPlaces {
    /// Create validated trap and hole collections without band-tail storage.
    pub fn new_standard(traps: Vec<Coord>, holes: Vec<Coord>) -> Result<Self, String> {
        for (index, trap) in traps.iter().enumerate() {
            trap.validate()
                .map_err(|error| format!("invalid trap {index}: {error}"))?;
        }
        for (index, hole) in holes.iter().enumerate() {
            hole.validate()
                .map_err(|error| format!("invalid hole {index}: {error}"))?;
        }

        Ok(Self::Standard { traps, holes })
    }

    /// Create collections from validated trap, hole, and bandtail coordinates.
    pub fn new_bandtail(
        traps: Vec<Coord>,
        holes: Vec<Coord>,
        bandtails: Vec<Coord>,
    ) -> Result<Self, String> {
        if bandtails.is_empty() {
            return Self::new_standard(traps, holes);
        }

        for (index, trap) in traps.iter().enumerate() {
            trap.validate()
                .map_err(|error| format!("invalid trap {index}: {error}"))?;
        }
        for (index, hole) in holes.iter().enumerate() {
            hole.validate()
                .map_err(|error| format!("invalid hole {index}: {error}"))?;
        }
        for (index, bandtail) in bandtails.iter().enumerate() {
            bandtail
                .validate()
                .map_err(|error| format!("invalid bandtail {index}: {error}"))?;
        }

        Ok(Self::WithBandtail {
            traps,
            holes,
            bandtails,
        })
    }

    /// Return all localised electron traps, regardless of variant.
    pub fn traps(&self) -> &[Coord] {
        match self {
            Self::Standard { traps, .. } | Self::WithBandtail { traps, .. } => traps,
        }
    }

    /// Return mutable trap storage independent of the outer layout variant.
    fn traps_mut(&mut self) -> &mut Vec<Coord> {
        match self {
            Self::Standard { traps, .. } | Self::WithBandtail { traps, .. } => traps,
        }
    }

    /// Return all hole recombination sites, regardless of variant.
    pub fn holes(&self) -> &[Coord] {
        match self {
            Self::Standard { holes, .. } | Self::WithBandtail { holes, .. } => holes,
        }
    }

    /// Return mutable hole storage independent of the outer layout variant.
    fn holes_mut(&mut self) -> &mut Vec<Coord> {
        match self {
            Self::Standard { holes, .. } | Self::WithBandtail { holes, .. } => holes,
        }
    }

    /// Return the band-tail sites when this is a band-tail simulation.
    pub fn bandtails(&self) -> Option<&[Coord]> {
        match self {
            Self::Standard { .. } => None,
            Self::WithBandtail { bandtails, .. } => Some(bandtails),
        }
    }

    /// Return mutable band-tail storage when the current layout provides it.
    fn bandtails_mut(&mut self) -> Option<&mut Vec<Coord>> {
        match self {
            Self::Standard { .. } => None,
            Self::WithBandtail { bandtails, .. } => Some(bandtails),
        }
    }

    /// Create empty site collections with capacity reserved for later assignment.
    pub fn with_capacity(
        traps: usize,
        holes: usize,
        bandtails: usize,
    ) -> Result<Self, String> {
        if bandtails == 0 {
            Ok(Self::Standard {
                traps: Self::reserved_vec(traps, "traps")?,
                holes: Self::reserved_vec(holes, "holes")?,
            })
        } else {
            Ok(Self::WithBandtail {
                traps: Self::reserved_vec(traps, "traps")?,
                holes: Self::reserved_vec(holes, "holes")?,
                bandtails: Self::reserved_vec(bandtails, "bandtails")?,
            })
        }
    }

    /// Create site collections populated with random coordinates.
    pub fn random_new<X: Numeric, Y: Numeric, Z: Numeric,>(
        t_no: usize,
        h_no: usize,
        b_no: usize,
        x: X,
        y: Y,
        z: Z,
        rng: & mut impl Rng,
    ) -> Result<Self, String> {
        let traps = (0..t_no)
            .map(|_| Coord::random_in(x, y, z, rng))
            .collect::<Result<Vec<_>, _>>()?;

        let holes = (0..h_no)
            .map(|_| Coord::random_in(x, y, z, rng))
            .collect::<Result<Vec<_>, _>>()?;
        if b_no > 0 {
            let bandtails = (0..b_no)
                .map(|_| Coord::random_in(x, y, z, rng))
                .collect::<Result<Vec<_>, _>>()?;
            ElectronPlaces::new_bandtail(traps, holes, bandtails)
        } else {
            ElectronPlaces::new_standard(traps, holes)
        }
    }

    /// Generate every site requested by a cube at an independent uniform position.
    ///
    /// The coordinates occupy the full rectangular volume. Correlations,
    /// exclusion radii, and crystallographic site constraints are not applied.
    pub fn random_from_cube(cube: &crate::crystal::Cube, rng: & mut impl Rng,) -> Result<Self, String> {
        let x = cube.boundary.x;
        let y = cube.boundary.y;
        let z = cube.boundary.z;
        let t_no = cube.trap_total;
        let h_no = cube.hole_total;
        let b_no = cube.bandtail_total;

        ElectronPlaces::random_new(t_no, h_no, b_no, x, y, z, rng)
    }

    /// Reserve exact fallible storage and attach a site-category name to errors.
    fn reserved_vec<T>(capacity: usize, name: &str) -> Result<Vec<T>, String> {
        let mut values = Vec::new();
        values
            .try_reserve_exact(capacity)
            .map_err(|error| format!("could not reserve capacity for {capacity} {name}: {error}"))?;
        Ok(values)
    }

    /// Replace all trap coordinates after validating them.
    pub fn set_traps(&mut self, traps: Vec<Coord>) -> Result<(), String> {
        for (index, trap) in traps.iter().enumerate() {
            trap.validate()
                .map_err(|error| format!("invalid trap {index}: {error}"))?;
        }
        *self.traps_mut() = traps;
        Ok(())
    }

    /// Replace all hole coordinates after validating them.
    pub fn set_holes(&mut self, holes: Vec<Coord>) -> Result<(), String> {
        for (index, hole) in holes.iter().enumerate() {
            hole.validate()
                .map_err(|error| format!("invalid hole {index}: {error}"))?;
        }
        *self.holes_mut() = holes;
        Ok(())
    }

    /// Replace all bandtail coordinates after validating them.
    pub fn set_bandtails(&mut self, bandtails: Vec<Coord>) -> Result<(), String> {
        for (index, bandtail) in bandtails.iter().enumerate() {
            bandtail.validate()
                .map_err(|error| format!("invalid bandtail {index}: {error}"))?;
        }
        match self {
            Self::Standard { traps, holes } if !bandtails.is_empty() => {
                *self = Self::WithBandtail {
                    traps: std::mem::take(traps),
                    holes: std::mem::take(holes),
                    bandtails,
                };
            }
            Self::Standard { .. } => {}
            Self::WithBandtail {
                traps,
                holes,
                bandtails: _,
            } if bandtails.is_empty() => {
                *self = Self::Standard {
                    traps: std::mem::take(traps),
                    holes: std::mem::take(holes),
                };
            }
            Self::WithBandtail {
                bandtails: current, ..
            } => *current = bandtails,
        }
        Ok(())
    }

    /// Push a pre-built trap into the collection.
    pub fn push_trap(&mut self, trap: Coord) -> Result<(), String> {
        trap.validate()?;
        self.traps_mut().push(trap);
        Ok(())
    }

    /// Create and append a trap at the given coordinate.
    pub fn push_trap_at<X: Numeric, Y: Numeric, Z: Numeric>(
        &mut self,
        x: X,
        y: Y,
        z: Z,
    ) -> Result<(), String> {
        self.push_trap(Coord::new(x, y, z)?)
    }

    /// Push a pre-built hole into the collection.
    pub fn push_hole(&mut self, hole: Coord) -> Result<(), String> {
        hole.validate()?;
        self.holes_mut().push(hole);
        Ok(())
    }

    /// Create and append a hole at the given coordinate.
    pub fn push_hole_at<X: Numeric, Y: Numeric, Z: Numeric>(
        &mut self,
        x: X,
        y: Y,
        z: Z,
    ) -> Result<(), String> {
        self.push_hole(Coord::new(x, y, z)?)
    }

    /// Push a pre-built bandtail into the collection.
    pub fn push_bandtail(&mut self, bandtail: Coord) -> Result<(), String> {
        bandtail.validate()?;
        match self.bandtails_mut() {
            Some(bandtails) => {
                bandtails.push(bandtail);
                Ok(())
            }
            None => self.set_bandtails(vec![bandtail]),
        }
    }

    /// Create and append a bandtail at the given coordinate.
    pub fn push_bandtail_at<X: Numeric, Y: Numeric, Z: Numeric>(
        &mut self,
        x: X,
        y: Y,
        z: Z,
    ) -> Result<(), String> {
        self.push_bandtail(Coord::new(x, y, z)?)
    }

    /// Calculate the Euclidean distance between two traps.
    ///
    /// This method does not apply a cube boundary. Use [`crate::crystal::Cube::distance`]
    /// with the returned trap coordinates when periodic wrapping is required.
    pub fn trap_trap_distance(&self, p1: usize, p2: usize) -> Float {
        self.traps()[p1].distance(&self.traps()[p2])
    }
    /// Replace every trap, hole, and bandtail coordinate with a random point.
    ///
    /// Site counts and the boundary remain unchanged. All new coordinates are
    /// generated from `0` through the corresponding boundary extent.
    pub fn randomise_positions<X: Numeric, Y: Numeric, Z: Numeric>(
        &mut self,
        x: &X,
        y: &Y,
        z: &Z,
        rng: & mut impl Rng,
    ) -> Result<(), String> {
        // Generate new random positions within the cube.
        for trap in self.traps_mut() {
            *trap = Coord::random_in(*x, *y, *z, rng)?;
        }

        for hole in self.holes_mut() {
            *hole = Coord::random_in(*x, *y, *z, rng)?;
        }

        if let Some(bandtails) = self.bandtails_mut() {
            for bandtail in bandtails {
                *bandtail = Coord::random_in(*x, *y, *z, rng)?;
            }
        }

        Ok(())
    }
}

/// Physical parameters assigned to one family of electron traps.
///
/// A trapped electron is represented by thermally coupled ground and excited
/// states. Each state has its own conduction-band activation energy,
/// delocalised attempt frequency, localised tunnelling prefactor, and spatial
/// decay constant. All energies are in electronvolts, frequencies and rates
/// are in inverse seconds, and inverse-length parameters must use units
/// reciprocal to the coordinate units.
#[derive(Debug, Clone, Copy)]
pub struct TrapParameters {
    /// Energy separation between the localised ground and excited states, in eV.
    pub excited_energy_gap: Float,
    /// Attempt frequency for thermal excitation from ground to excited state, in s⁻¹.
    pub s_frequency_e: Float,
    /// Relaxation frequency from excited to ground state, in s⁻¹.
    pub s_frequency_g: Float,
    /// Ground-state activation energy for release to the conduction band, in eV.
    pub e_cb_ground: Float,
    /// Excited-state activation energy for release to the conduction band, in eV.
    pub e_cb_excited: Float,
    /// Ground-state delocalised-release attempt frequency, in s⁻¹.
    pub de_frequency_ground: Float,
    /// Excited-state delocalised-release attempt frequency, in s⁻¹.
    pub de_frequency_excited: Float,
    /// Ground-state localised-tunnelling prefactor, in s⁻¹.
    pub lo_frequency_ground: Float,
    /// Excited-state localised-tunnelling prefactor, in s⁻¹.
    pub lo_frequency_excited: Float,
    /// Ground-state tunnelling decay constant in inverse coordinate units.
    pub alpha_ground: Float,
    /// Excited-state tunnelling decay constant in inverse coordinate units.
    pub alpha_excited: Float,
    /// Delocalised destination length scale in the same units as coordinates.
    pub delocalised_mu: Float,
    /// Factor multiplying the delocalised retrapping mean waiting time.
    ///
    /// With the current reciprocal-rate model, larger positive values make
    /// retrapping slower relative to the unit recombination factor; zero
    /// disables conduction-band retrapping destinations.
    pub retrap_ratio: Float,
}
impl TrapParameters {
    /// Construct a complete trap-family parameter record without validation.
    ///
    /// The caller must keep energies, frequencies, and length scales in the
    /// units described by [`TrapParameters`] and ensure they are physically valid.
    pub fn new(
        excited_energy_gap: Float,
        s_frequency_e: Float,
        s_frequency_g: Float,
        e_cb_ground: Float,
        e_cb_excited: Float,
        de_frequency_ground: Float,
        de_frequency_excited: Float,
        lo_frequency_ground: Float,
        lo_frequency_excited: Float,
        alpha_ground: Float,
        alpha_excited: Float,
        delocalised_mu: Float,
        retrap_ratio: Float,
    ) -> Self {
        Self {
            excited_energy_gap,
            s_frequency_e,
            s_frequency_g,
            e_cb_ground,
            e_cb_excited,
            de_frequency_ground,
            de_frequency_excited,
            lo_frequency_ground,
            lo_frequency_excited,
            alpha_ground,
            alpha_excited,
            delocalised_mu,
            retrap_ratio,
        }
    }
}

/// Maps spatial trap identifiers to their physical parameter records.
///
/// Uniform layouts avoid duplicating a record for homogeneous crystals,
/// while direct and indexed layouts support heterogeneous trap populations.
pub enum TrapParameterLayout {
    /// Every trap uses exactly the same parameters.
    Uniform(TrapParameters),

    /// Every trap has its own parameters.
    ///
    /// `parameters[trap_id]`
    Direct(Box<[TrapParameters]>),

    /// Traps reference a table of shared parameter records.
    ///
    /// Useful for families and mixtures of shared/individual parameters.
    Indexed {
        /// Table of distinct parameter records shared between traps.
        records: Box<[TrapParameters]>,
        /// Parameter-record identifier assigned to each trap identifier.
        by_trap: Box<[PlaceId]>,
    },
}

impl TrapParameterLayout {
    /// Build a homogeneous layout from the individual model parameters.
    ///
    /// The excited-state conduction-band barrier is derived as
    /// `e_cb - excited_energy_gap`; every trap subsequently borrows the same
    /// resulting [`TrapParameters`] record.
    pub fn new_uniform(excited_energy_gap: Float, s_frequency_e: Float, 
                       s_frequency_g: Float, e_cb: Float, 
                       s_gs: Float, s_es: Float, 
                       b_gs: Float, b_es: Float, 
                       alpha_gs: Float, alpha_es: Float, 
                       mu: Float, retrap_ratio: Float) -> Self {

        let parameters = TrapParameters::new(excited_energy_gap,
                                                             s_frequency_e,
                                                             s_frequency_g,
                                                             e_cb,
                                                             e_cb - excited_energy_gap,
                                                             s_gs,
                                                             s_es,
                                                             b_gs,
                                                             b_es,
                                                             alpha_gs,
                                                             alpha_es,
                                                             mu,
                                                             retrap_ratio,);
        TrapParameterLayout::uniform(parameters)
    }

    /// Assign one parameter record to every trap identifier.
    pub fn uniform(parameters: TrapParameters) -> Self {
        Self::Uniform(parameters)
    }

    /// Assign one parameter record per trap in identifier order.
    ///
    /// Returns an error unless `parameters.len()` equals `trap_count`.
    pub fn direct(parameters: Vec<TrapParameters>, trap_count: usize) -> Result<Self, String> {
        if parameters.len() != trap_count {
            return Err(format!(
                "expected {trap_count} trap parameter records, found {}",
                parameters.len(),
            ));
        }

        Ok(Self::Direct(parameters.into_boxed_slice()))
    }

    /// Assign traps to a shared table of parameter records.
    ///
    /// `assignments[trap.index()]` identifies the entry in `records`. The
    /// constructor validates both the trap count and every table reference.
    pub fn indexed(
        records: Vec<TrapParameters>,
        assignments: Vec<PlaceId>,
        trap_count: usize,
    ) -> Result<Self, String> {
        if records.is_empty() {
            return Err("at least one parameter record is required".to_string());
        }

        if assignments.len() != trap_count {
            return Err(format!(
                "expected {trap_count} parameter assignments, found {}",
                assignments.len(),
            ));
        }

        for (trap_index, parameter_id) in assignments.iter().enumerate() {
            if parameter_id.index() >= records.len() {
                return Err(format!(
                    "trap {trap_index} references missing parameter {}",
                    parameter_id.index(),
                ));
            }
        }

        Ok(Self::Indexed {
            records: records.into_boxed_slice(),
            by_trap: assignments.into_boxed_slice(),
        })
    }

    /// Return the physical parameters assigned to `trap`.
    ///
    /// Layout constructors guarantee valid indexing; passing a trap ID from a
    /// different crystal can therefore panic.
    pub fn get(&self, trap: PlaceId) -> &TrapParameters {
        match self {
            Self::Uniform(parameters) => parameters,

            Self::Direct(parameters) => &parameters[trap.index()],

            Self::Indexed { records, by_trap } => &records[by_trap[trap.index()].index()],
        }
    }
    /// Returns the DelocalisedTransitionInputs for all the available traps
    pub fn get_delocalised<'a>(
        &self,
        ground_weights: &'a [Float],
        excited_weights: &'a [Float],
        traps: &PlaceAvailability,
        temperature: Float,
        ) -> Result<
        DelocalisedTransitionInputs<Vec<Float>, Vec<Float>, &'a [Float]>,
        String,> 
    {
        let available = traps.available();

        if ground_weights.len() != available.len()
            || excited_weights.len() != available.len()
        {
            return Err(format!(
                "expected {} weights, found {} ground and {} excited",
                available.len(),
                ground_weights.len(),
                excited_weights.len(),
            ));
        }

        let e_cb_ground =  match self {
            Self::Uniform(parameters) => {
                vec![parameters.e_cb_ground]
            }, 
            Self::Direct(parameters) => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    parameters[index]
                    .e_cb_ground).collect::<Vec<Float>>()
            },
            Self::Indexed { records, by_trap } => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    records[by_trap[index].index()]
                    .e_cb_ground).collect::<Vec<Float>>()
            }
        };
      
        let frequency_ground =  match self {
            Self::Uniform(parameters) => {
                vec![parameters.de_frequency_ground]
            }, 
            Self::Direct(parameters) => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    parameters[index]
                    .de_frequency_ground).collect::<Vec<Float>>()
            },
            Self::Indexed { records, by_trap } => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    records[by_trap[index].index()]
                    .de_frequency_ground).collect::<Vec<Float>>()
            }
        };
        let e_cb_excited = match self {
            Self::Uniform(parameters) => {
                vec![parameters.e_cb_excited]
            }, 
            Self::Direct(parameters) => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    parameters[index]
                    .e_cb_excited).collect::<Vec<Float>>()
            },
            Self::Indexed { records, by_trap } => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    records[by_trap[index].index()]
                    .e_cb_excited).collect::<Vec<Float>>()
            }
        };
        
        let frequency_excited = match self {
            Self::Uniform(parameters) => {
                vec![parameters.de_frequency_excited]
            }, 
            Self::Direct(parameters) => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    parameters[index]
                    .de_frequency_excited).collect::<Vec<Float>>()
            },
            Self::Indexed { records, by_trap } => {
                traps.available_indices_vec()
                    .into_iter()
                    .map(|index| 
                    records[by_trap[index].index()]
                    .de_frequency_excited).collect::<Vec<Float>>()
            }
        };
        
        Ok(DelocalisedTransitionInputs {
            e_cb_ground,
            frequency_ground,
            e_cb_excited,
            frequency_excited,
            temperature,
            ground_weight: ground_weights,
            excited_weight: excited_weights,
        })
    }

       /// Returns the DelocalisedTransitionInputs for all the available traps
    pub fn get_localised_recombination<'a>(
        &self,
        ground_weight: &'a Float,
        excited_weight: &'a Float,
        distance: &'a [Float],
        trap: &PlaceId,
        temperature: Float,
    ) -> Result<
        LocalisedTransitionInputs<Float, Float, &'a Float, &'a [Float]>,
        String,> 
    {
        Ok( match self {
            Self::Uniform(parameters) => {
                LocalisedTransitionInputs {
                    alpha_ground: parameters.alpha_ground,
                    frequency_ground: parameters.lo_frequency_ground,
                    alpha_excited: parameters.alpha_excited,
                    frequency_excited: parameters.lo_frequency_excited,
                    ground_weight,
                    excited_weight,
                    distance,
                }
            },
             Self::Direct(parameters) => { LocalisedTransitionInputs {
                    alpha_ground: parameters[trap.index()].alpha_ground,
                    frequency_ground: parameters[trap.index()].lo_frequency_ground,
                    alpha_excited: parameters[trap.index()].alpha_excited,
                    frequency_excited: parameters[trap.index()].lo_frequency_excited,
                    ground_weight,
                    excited_weight,
                    distance,
                }
             },

            Self::Indexed { records, by_trap } => {LocalisedTransitionInputs {
            alpha_ground: records[by_trap[trap.index()].index()].alpha_ground,
                    frequency_ground: records[by_trap[trap.index()].index()].lo_frequency_ground,
                    alpha_excited: records[by_trap[trap.index()].index()].alpha_excited,
                    frequency_excited: records[by_trap[trap.index()].index()].lo_frequency_excited,
                    ground_weight,
                    excited_weight,
                    distance,
                }
             },
        
        })
    }

    pub fn get_localised_retrapping<'a>(
        &self,
        ground_weight: &'a Float,
        excited_weight: &'a Float,
        distance: &'a [Float],
        trap: &PlaceId,
        temperature: Float,
    ) -> Result<
        LocalisedTransitionInputs<Float, Float, &'a Float, &'a [Float]>,
        String,> 
    {
        Ok( match self {
            Self::Uniform(parameters) => {
                LocalisedTransitionInputs {
                    alpha_ground: parameters.alpha_ground,
                    frequency_ground: parameters.lo_frequency_ground,
                    alpha_excited: parameters.alpha_excited,
                    frequency_excited: parameters.lo_frequency_excited,
                    ground_weight,
                    excited_weight,
                    distance,
                }
            },
             Self::Direct(parameters) => { LocalisedTransitionInputs {
                    alpha_ground: parameters[trap.index()].alpha_ground,
                    frequency_ground: parameters[trap.index()].lo_frequency_ground,
                    alpha_excited: parameters[trap.index()].alpha_excited,
                    frequency_excited: parameters[trap.index()].lo_frequency_excited,
                    ground_weight,
                    excited_weight,
                    distance,
                }
             },

            Self::Indexed { records, by_trap } => {LocalisedTransitionInputs {
            alpha_ground: records[by_trap[trap.index()].index()].alpha_ground,
                    frequency_ground: records[by_trap[trap.index()].index()].lo_frequency_ground,
                    alpha_excited: records[by_trap[trap.index()].index()].alpha_excited,
                    frequency_excited: records[by_trap[trap.index()].index()].lo_frequency_excited,
                    ground_weight,
                    excited_weight,
                    distance,
                }
             },
        
        })

       
        
       
    }

}


#[cfg(test)]
mod tests {
    use crate::random;

use super::*;
    fn print_coordinates(label: &str, coordinates: &[Coord]) {
        println!("{label}:");
        for (index, coordinate) in coordinates.iter().enumerate() {
            println!(
                "  {index}: x = {}, y = {}, z = {}",
                coordinate.x, coordinate.y, coordinate.z,
            );
        }
    }

    fn every_position_changed(original: &[Coord], updated: &[Coord]) -> bool {
        original.len() == updated.len()
            && original.iter().zip(updated).all(|(original, updated)| {
                original.x != updated.x
                    || original.y != updated.y
                    || original.z != updated.z
            })
    }
     #[test]
    fn generation_rejects_invalid_coordinates() {
        let mut rng = random::get_std_rng_for_rep(0);
        assert!(Coord::new(-1.0, 0.0, 0.0).is_err());
        assert!(Coord::random_in(10.0, -1.0, 10.0, &mut rng).is_err());
    }

    #[test]
    fn distance_euclidean() {
        let p1 = Coord::new(0.0, 0.0, 0.0).unwrap();
        let p2 = Coord::new(3.0, 4.0, 0.0).unwrap();
        assert_eq!(p1.distance(&p2), 5.0);
    }


    #[test]
    fn randomising_cube_positions_changes_traps_holes_and_bandtails() {
        let mut rng = random::get_std_rng_for_rep(0);
        let mut places = ElectronPlaces::random_new(3, 3, 3, 10.0, 10.0, 10.0, &mut rng).unwrap();
        let cube = crate::crystal::Cube::new(10.0,10.0,10.0,3,3,3,true).unwrap();
        let original_places = places.clone();

        println!("Original coordinates");
        print_coordinates("Traps", original_places.traps());
        print_coordinates("Holes", original_places.holes());
        print_coordinates("Bandtails", original_places.bandtails().unwrap());

        places.randomise_positions(&cube.boundary.x, &cube.boundary.y, &cube.boundary.z, &mut rng).unwrap();

        println!("Randomised coordinates");
        print_coordinates("Traps", places.traps());
        print_coordinates("Holes", places.holes());
        print_coordinates("Bandtails", places.bandtails().unwrap());

        assert!(every_position_changed(
            original_places.traps(),
            places.traps(),
        ));
        assert!(every_position_changed(
            original_places.holes(),
            places.holes(),
        ));
        assert!(every_position_changed(
            original_places.bandtails().unwrap(),
            places.bandtails().unwrap(),
        ));

        assert!(places.traps().iter().all(|trap| cube.contains(trap)));
        assert!(places.holes().iter().all(|hole| cube.contains(hole)));
        assert!(places
            .bandtails()
            .unwrap()
            .iter()
            .all(|bandtail| cube.contains(bandtail)));
    }

     #[test]
    fn site_generation_propagates_coordinate_errors() {
        let invalid = Coord {
            x: Float::NAN,
            y: 0.0,
            z: 0.0,
        };

        assert!(invalid.validate().is_err());
    }

    #[test]
    fn electron_place_mutations_validate_before_changing_collections() {
        let mut places = ElectronPlaces::with_capacity(2, 2, 2).unwrap();
        places.push_trap_at(1.0, 2.0, 3.0).unwrap();
        places.push_hole_at(2.0, 3.0, 4.0).unwrap();
        places.push_bandtail_at(3.0, 4.0, 5.0).unwrap();

        assert_eq!(places.traps().len(), 1);
        assert_eq!(places.holes().len(), 1);
        assert_eq!(places.bandtails().unwrap().len(), 1);

        places.set_traps(vec![Coord::new(4.0, 5.0, 6.0).unwrap()]).unwrap();
        places.set_holes(vec![Coord::new(5.0, 6.0, 7.0).unwrap()]).unwrap();
        places
            .set_bandtails(vec![Coord::new(6.0, 7.0, 8.0).unwrap()])
            .unwrap();

        let invalid = Coord {
            x: -1.0,
            y: 0.0,
            z: 0.0,
        };
        let error = places.set_traps(vec![invalid]).unwrap_err();
        assert!(error.contains("invalid trap 0"));
        assert_eq!(places.traps().len(), 1);
        assert_eq!(places.traps()[0].x, 4.0);
    }

    #[test]
    fn constructors_reserve_and_populate_requested_counts_and_report_allocation_errors() {
        let empty = ElectronPlaces::with_capacity(3, 6, 9)
                                    .expect("small capacities should be reservable");
        match empty {
            ElectronPlaces::WithBandtail {
                traps,
                holes,
                bandtails,
            } => {
                assert_eq!(traps.capacity(), 3);
                assert_eq!(holes.capacity(), 6);
                assert_eq!(bandtails.capacity(), 9);
            }
            ElectronPlaces::Standard { .. } => panic!("expected band-tail storage"),
        }
        let mut rng = random::get_std_rng_for_rep(0);
        let populated = ElectronPlaces::random_new(3, 6, 9, 2.0, 1.0, 1.0, &mut rng)
            .expect("small random cube should be generated");
        assert_eq!(populated.traps().len(), 3);
        assert_eq!(populated.holes().len(), 6);
        assert_eq!(populated.bandtails().unwrap().len(), 9);

        let error = ElectronPlaces::with_capacity(2, usize::MAX, 0).unwrap_err();
        assert!(
            error.starts_with(&format!(
                "could not reserve capacity for {} holes:",
                usize::MAX
            )),
            "unexpected allocation error: {error}"
        );
    }

    #[test]
    fn standard_places_add_and_remove_bandtail_storage() {
        let mut places = ElectronPlaces::with_capacity(2, 2, 0).unwrap();
        assert!(matches!(places, ElectronPlaces::Standard { .. }));
        assert!(places.bandtails().is_none());

        places.push_bandtail_at(1.0, 2.0, 3.0).unwrap();
        assert!(matches!(places, ElectronPlaces::WithBandtail { .. }));
        assert_eq!(places.bandtails().unwrap().len(), 1);

        places.set_bandtails(Vec::new()).unwrap();
        assert!(matches!(places, ElectronPlaces::Standard { .. }));
        assert!(places.bandtails().is_none());
    }
}
