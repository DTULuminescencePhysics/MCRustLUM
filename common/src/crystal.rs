// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Spatial representation of the simulation crystal.
//!
//! A [`crate::crystal::Cube`] combines its boundary condition with the
//! positions of electron traps, hole recombination sites, and bandtail states.
//! Constructors validate geometry before allocating or randomly generating
//! site coordinates.
use crate::numeric::{Float, Numeric};

/// Distance behaviour at the faces of a simulation volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    /// Opposite faces are connected when calculating distances.
    Periodic,
    /// Coordinates do not wrap across opposite faces.
    Padded,
}

/// Positive x, y, and z extents together with a boundary condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boundary {
    /// Maximum x coordinate.
    pub x: Float,
    /// Maximum y coordinate.
    pub y: Float,
    /// Maximum z coordinate.
    pub z: Float,
    /// Distance behavior at the cube faces.
    pub kind: BoundaryCondition,
}

impl Boundary {
    /// Create a Periodic or Padded boundary condition.
    pub fn new<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        periodic: bool,
    ) -> Result<Self, String> {
        let x = x.to_float();
        let y = y.to_float();
        let z = z.to_float();

        Self::validate_dimension("x", x)?;
        Self::validate_dimension("y", y)?;
        Self::validate_dimension("z", z)?;

        let kind = if periodic {
            BoundaryCondition::Periodic
        } else {
            BoundaryCondition::Padded
        };

        Ok(Self { x, y, z, kind })
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
    /// Check that every extent is finite and greater than zero.
    pub fn validate(&self) -> Result<(), String> {
        Self::validate_dimension("x", self.x)?;
        Self::validate_dimension("y", self.y)?;
        Self::validate_dimension("z", self.z)
    }

    /// Calculate the ordinary Euclidean distance between two points.
    pub fn padded_distance(p1: &Coord, p2: &Coord) -> Float {
        ((p1.x - p2.x).powi(2) + (p1.y - p2.y).powi(2) + (p1.z - p2.z).powi(2)).sqrt()
    }

    /// Calculate the shortest distance after wrapping across opposite faces.
    pub fn periodic_distance(&self, p1: &Coord, p2: &Coord) -> Float {
        let dx = if (p1.x - p2.x).abs() > self.x / 2.0 {
            self.x - (p1.x - p2.x).abs()
        } else {
            (p1.x - p2.x).abs()
        };
        let dy = if (p1.y - p2.y).abs() > self.y / 2.0 {
            self.y - (p1.y - p2.y).abs()
        } else {
            (p1.y - p2.y).abs()
        };
        let dz = if (p1.z - p2.z).abs() > self.z / 2.0 {
            self.z - (p1.z - p2.z).abs()
        } else {
            (p1.z - p2.z).abs()
        };
        (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt()
    }

    /// Calculate boundary-aware distance between two coordinates.
    pub fn distance(&self, p1: &Coord, p2: &Coord) -> Float {
        match self.kind {
            BoundaryCondition::Padded => Boundary::padded_distance(p1, p2),
            BoundaryCondition::Periodic => self.periodic_distance(p1, p2),
        }
    }
    /// Return whether a coordinate lies inside or on the boundary.
    pub fn contains(&self, point: &Coord) -> bool {
        point.x >= 0.0
            && point.x <= self.x
            && point.y >= 0.0
            && point.y <= self.y
            && point.z >= 0.0
            && point.z <= self.z
    }
    /// Return `x * y * z` in the cube's length units cubed.
    pub fn volume(&self) -> Float {
        self.x * self.y * self.z
    }
    /// Convert a volumetric site density to a whole number of sites.
    ///
    /// The positive result is truncated toward zero. An error is returned for
    /// a negative or non-finite density and when the result exceeds `usize`.
    pub fn density_to_number(&self, density: Float) -> Result<usize, String> {
        self.validate()?;

        if !density.is_finite() {
            return Err(format!("density must be finite, got {density}"));
        }
        if density < 0.0 {
            return Err(format!("density cannot be negative, got {density}"));
        }

        let count = density * self.volume();
        if !count.is_finite() || count > usize::MAX as Float {
            return Err(format!("density produces an invalid trap count: {count}"));
        }

        Ok(count as usize)
    }
}

/// A three-dimensional site position.
///
/// Coordinates can be supplied directly or generated within positive x, y,
/// and z limits.
#[derive(Debug, Clone)]
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
    pub fn random_in<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
    ) -> Result<Self, String> {
        let x = x.to_float();
        let y = y.to_float();
        let z = z.to_float();

        Boundary::validate_dimension("x", x)?;
        Boundary::validate_dimension("y", y)?;
        Boundary::validate_dimension("z", z)?;

        let mut rng = rand::thread_rng();
        Ok(Coord {
            x: Float::random_in(x, &mut rng),
            y: Float::random_in(y, &mut rng),
            z: Float::random_in(z, &mut rng),
        })
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
pub struct ElectronPlaces {
    /// Positions of localised electron traps.
    pub traps: Vec<Coord>,
    /// Positions of hole recombination sites.
    pub holes: Vec<Coord>,
    /// Positions of shallow bandtail states.
    pub bandtails: Vec<Coord>,
}

impl ElectronPlaces {
    /// Create collections from validated trap, hole, and bandtail coordinates.
    pub fn new(
        traps: Vec<Coord>,
        holes: Vec<Coord>,
        bandtails: Vec<Coord>,
    ) -> Result<Self, String> {
        for (index, trap) in traps.iter().enumerate() {
            trap.validate()
                .map_err(|error| format!("invalid trap {index}: {error}"))?;
        }
        for (index, hole) in holes.iter().enumerate() {
            hole.validate()
                .map_err(|error| format!("invalid hole {index}: {error}"))?;
        }
        for (index, bandtail) in bandtails.iter().enumerate() {
            bandtail.validate()
                    .map_err(|error| format!("invalid bandtail {index}: {error}"))?;
        }

        Ok(Self {
            traps,
            holes,
            bandtails,
        })
    }

    /// Create empty site collections with capacity reserved for later assignment.
    pub fn with_capacity(
        traps: usize,
        holes: usize,
        bandtails: usize,
    ) -> Result<Self, String> {
        Ok(Self {
            traps: Self::reserved_vec(traps, "traps")?,
            holes: Self::reserved_vec(holes, "holes")?,
            bandtails: Self::reserved_vec(bandtails, "bandtails")?,
        })
    }

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
        self.traps = traps;
        Ok(())
    }

    /// Replace all hole coordinates after validating them.
    pub fn set_holes(&mut self, holes: Vec<Coord>) -> Result<(), String> {
        for (index, hole) in holes.iter().enumerate() {
            hole.validate()
                .map_err(|error| format!("invalid hole {index}: {error}"))?;
        }
        self.holes = holes;
        Ok(())
    }

    /// Replace all bandtail coordinates after validating them.
    pub fn set_bandtails(&mut self, bandtails: Vec<Coord>) -> Result<(), String> {
        for (index, bandtail) in bandtails.iter().enumerate() {
            bandtail.validate()
                .map_err(|error| format!("invalid bandtail {index}: {error}"))?;
        }
        self.bandtails = bandtails;
        Ok(())
    }

    /// Push a pre-built trap into the collection.
    pub fn push_trap(&mut self, trap: Coord) -> Result<(), String> {
        trap.validate()?;
        self.traps.push(trap);
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
        self.holes.push(hole);
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
        self.bandtails.push(bandtail);
        Ok(())
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
    /// This method does not apply a cube boundary. Use
    /// [`Cube::trap_trap_distance`] when periodic wrapping is required.
    pub fn trap_trap_distance(&self, p1: usize, p2: usize) -> Float {
        self.traps[p1].distance(&self.traps[p2])
    }

    /// Create site collections populated with random coordinates.
    pub fn random_new<X: Numeric, Y: Numeric, Z: Numeric>(
        t_no: usize,
        h_no: usize,
        b_no: usize,
        x: X,
        y: Y,
        z: Z,
    ) -> Result<Self, String> {
        let traps = (0..t_no)
            .map(|_| Coord::random_in(x, y, z))
            .collect::<Result<Vec<_>, _>>()?;

        let holes = (0..h_no)
            .map(|_| Coord::random_in(x, y, z))
            .collect::<Result<Vec<_>, _>>()?;

        let bandtails = (0..b_no)
            .map(|_| Coord::random_in(x, y, z))
            .collect::<Result<Vec<_>, _>>()?;

        ElectronPlaces::new(traps, holes, bandtails)
    }
}

/// A simulation volume containing all electron-site coordinates.
///
/// Construction validates dimensions and propagates random-generation errors.
/// Use a mutable cube to regenerate every stored position between Monte Carlo
/// iterations.
///
/// # Example
///
/// ```
/// # fn main() -> Result<(), String> {
/// use common::crystal::Cube;
///
/// let mut cube = Cube::new_random(10.0, 10.0, 10.0, 3, 3, 3, true)?;
/// assert_eq!(cube.places.traps.len(), 3);
/// cube.randomise_positions()?;
/// assert!(cube.places.traps.iter().all(|point| cube.contains(point)));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct Cube {
    /// Trap, hole, and bandtail coordinates in this realization.
    pub places: ElectronPlaces,
    /// Dimensions and boundary condition used by distance calculations.
    pub boundary: Boundary,
}

impl Cube {
    /// Create a validated cube with empty, preallocated site vectors.
    pub fn new_empty<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        t_no: usize,
        h_no: usize,
        b_no: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let boundary = Boundary::new(x, y, z, periodic)?;
        let places = ElectronPlaces::with_capacity(t_no, h_no, b_no)?;
        Ok(Self { places, boundary })
    }
    /// Create an empty cube whose capacities are derived from trap density.
    ///
    /// `h_no` and `b_no` are counts per trap, not absolute counts.
    pub fn new_empty_from_density<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        density: Float,
        h_no: usize,
        b_no: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let boundary = Boundary::new(x, y, z, periodic)?;
        let t_no = boundary.density_to_number(density)?;
        let h_no = h_no
            .checked_mul(t_no)
            .ok_or_else(|| "hole count overflowed usize".to_string())?;
        let b_no = b_no
            .checked_mul(t_no)
            .ok_or_else(|| "bandtail count overflowed usize".to_string())?;

        let places = ElectronPlaces::with_capacity(t_no, h_no, b_no)?;
        Ok(Self { places, boundary })
    }
    /// Create a cube with explicit numbers of randomly positioned sites.
    pub fn new_random<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        t_no: usize,
        h_no: usize,
        b_no: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let boundary = Boundary::new(x, y, z, periodic)?;
        let places = ElectronPlaces::random_new(
            t_no,
            h_no,
            b_no,
            boundary.x,
            boundary.y,
            boundary.z,
        )?;
        Ok(Self { places, boundary })
    }
    /// Create a cube with random sites derived from density and per-trap ratios.
    ///
    /// `h_no` and `b_no` are multiplied by the calculated trap count. Errors
    /// report invalid dimensions or density, arithmetic overflow, and failed
    /// allocation or coordinate generation.
    pub fn new_random_from_density<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        density: Float,
        h_no: usize,
        b_no: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let boundary = Boundary::new(x, y, z, periodic)?;
        let t_no = boundary.density_to_number(density)?;
        let h_no = h_no
            .checked_mul(t_no)
            .ok_or_else(|| "hole count overflowed usize".to_string())?;
        let b_no = b_no
            .checked_mul(t_no)
            .ok_or_else(|| "bandtail count overflowed usize".to_string())?;
        let places = ElectronPlaces::random_new(
            t_no,
            h_no,
            b_no,
            boundary.x,
            boundary.y,
            boundary.z,
        )?;
        Ok(Self { places, boundary })
    }

    /// Return whether `point` lies within this cube.
    pub fn contains(&self, point: &Coord) -> bool {
        self.boundary.contains(point)
    }

    /// Calculate the distance between two points using the cube's boundary condition.
    pub fn distance(&self, p1: &Coord, p2: &Coord) -> Float {
        self.boundary.distance(p1, p2)
    }

    /// Return the boundary-aware distance between two indexed traps.
    ///
    /// # Panics
    ///
    /// Panics if either index is outside the trap vector.
    pub fn trap_trap_distance(&self, p1: usize, p2: usize) -> Float {
        self.boundary.distance(
            &self.places.traps[p1],
            &self.places.traps[p2],
        )
    }

    /// Generate one random point inside the cube.
    pub fn random_point(&self) -> Result<Coord, String> {
        Coord::random_in(self.boundary.x, self.boundary.y, self.boundary.z)
    }
    /// Replace every trap, hole, and bandtail coordinate with a random point.
    ///
    /// Site counts and the boundary remain unchanged. All new coordinates are
    /// generated from `0` through the corresponding boundary extent.
    pub fn randomise_positions(&mut self) -> Result<(), String> {
        self.boundary.validate()?;

        let x = self.boundary.x;
        let y = self.boundary.y;
        let z = self.boundary.z;

        // Generate new random positions within the cube.
        for trap in &mut self.places.traps {
            *trap =  Coord::random_in(x, y, z)?;
        }

        for hole in &mut self.places.holes {
            *hole =  Coord::random_in(x, y, z)?;
        }

        for bandtail in &mut self.places.bandtails {
            *bandtail = Coord::random_in(x, y, z)?;
        }

        Ok(())
    }

}

#[cfg(test)]
mod tests {
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
    fn distance_euclidean() {
        let p1 = Coord::new(0.0, 0.0, 0.0).unwrap();
        let p2 = Coord::new(3.0, 4.0, 0.0).unwrap();
        assert_eq!(p1.distance(&p2), 5.0);
    }

    #[test]
    fn distance_padded_boundary() {
        let p1 = Coord::new(1.0, 1.0, 1.0).unwrap();
        let p2 = Coord::new(4.0, 1.0, 1.0).unwrap();
        let boundary = Boundary::new(10.0, 10.0, 10.0, false).unwrap();
        let distance = boundary.distance(&p1, &p2);
        assert_eq!(distance, 3.0);
    }

    #[test]
    fn distance_periodic_boundary_wraps() {
        let p1 = Coord::new(0.5, 5.0, 5.0).unwrap();
        let p2 = Coord::new(8.5, 5.0, 5.0).unwrap();
        let boundary = Boundary::new(10.0, 10.0, 10.0, true).unwrap();
        // Without wrapping: distance = 9.0
        // With wrapping: distance should be 2.0 (shorter path wraps around)
        let distance = boundary.distance(&p1, &p2);
        assert_eq!(distance, 2.0);
    }

    #[test]
    fn cube_distance_with_boundary() {
        let mut cube = Cube::new_random(10.0, 10.0, 10.0, 0, 0, 0, true).unwrap();
        cube.places.push_trap_at(1.0, 5.0, 5.0).unwrap();
        cube.places.push_trap_at(9.0, 5.0, 5.0).unwrap();
        // Periodic: shortest distance is 2.0
        assert_eq!(cube.trap_trap_distance(0, 1), 2.0);
    }

    #[test]
    fn randomising_cube_positions_changes_traps_holes_and_bandtails() {
        let mut cube = Cube::new_random(10.0, 10.0, 10.0, 3, 3, 3, true).unwrap();
        let original_places = cube.places.clone();

        println!("Original coordinates");
        print_coordinates("Traps", &original_places.traps);
        print_coordinates("Holes", &original_places.holes);
        print_coordinates("Bandtails", &original_places.bandtails);

        cube.randomise_positions().unwrap();

        println!("Randomised coordinates");
        print_coordinates("Traps", &cube.places.traps);
        print_coordinates("Holes", &cube.places.holes);
        print_coordinates("Bandtails", &cube.places.bandtails);

        assert!(every_position_changed(
            &original_places.traps,
            &cube.places.traps,
        ));
        assert!(every_position_changed(
            &original_places.holes,
            &cube.places.holes,
        ));
        assert!(every_position_changed(
            &original_places.bandtails,
            &cube.places.bandtails,
        ));

        assert!(cube.places.traps.iter().all(|trap| cube.contains(trap)));
        assert!(cube.places.holes.iter().all(|hole| cube.contains(hole)));
        assert!(cube
            .places
            .bandtails
            .iter()
            .all(|bandtail| cube.contains(bandtail)));
    }

    #[test]
    fn coord_with_integers() {
        let boundary = Boundary::new(10i32, 10i32, 10i32, false).unwrap();
        let p1 = Coord::new(0i32, 0i32, 0i32).unwrap();
        let p2 = Coord::new(3i32, 4i32, 0i32).unwrap();
        boundary.distance(&p1, &p2);
        assert_eq!(p1.distance(&p2), 5.0);
    }

    #[test]
    fn generation_rejects_invalid_boundaries_coordinates_and_density() {
        assert!(Boundary::new(0.0, 10.0, 10.0, false).is_err());
        assert!(Boundary::new(Float::NAN, 10.0, 10.0, false).is_err());
        assert!(Coord::new(-1.0, 0.0, 0.0).is_err());
        assert!(Coord::random_in(10.0, -1.0, 10.0).is_err());
        assert!(Cube::new_random_from_density(
            10.0,
            10.0,
            10.0,
            -1.0,
            1,
            1,
            false,
        )
        .is_err());
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
    fn boundary_volume_density_and_containment_include_the_faces() {
        let boundary = Boundary::new(2.0, 3.0, 4.0, false).unwrap();

        assert_eq!(boundary.volume(), 24.0);
        assert_eq!(boundary.density_to_number(0.5).unwrap(), 12);
        assert!(boundary.contains(&Coord::new(0.0, 0.0, 0.0).unwrap()));
        assert!(boundary.contains(&Coord::new(2.0, 3.0, 4.0).unwrap()));
        assert!(!boundary.contains(&Coord::new(2.1, 3.0, 4.0).unwrap()));

        let invalid_boundary = Boundary {
            x: 2.0,
            y: Float::INFINITY,
            z: 4.0,
            kind: BoundaryCondition::Padded,
        };
        assert!(invalid_boundary.validate().is_err());
        assert!(boundary.density_to_number(Float::MAX).is_err());
    }

    #[test]
    fn electron_place_mutations_validate_before_changing_collections() {
        let mut places = ElectronPlaces::with_capacity(2, 2, 2).unwrap();
        places.push_trap_at(1.0, 2.0, 3.0).unwrap();
        places.push_hole_at(2.0, 3.0, 4.0).unwrap();
        places.push_bandtail_at(3.0, 4.0, 5.0).unwrap();

        assert_eq!(places.traps.len(), 1);
        assert_eq!(places.holes.len(), 1);
        assert_eq!(places.bandtails.len(), 1);

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
        assert_eq!(places.traps.len(), 1);
        assert_eq!(places.traps[0].x, 4.0);
    }

    #[test]
    fn density_constructors_apply_per_trap_ratios_and_report_overflow() {
        let empty = Cube::new_empty_from_density(2.0, 1.0, 1.0, 1.5, 2, 3, false)
            .expect("small capacities should be reservable");
        assert_eq!(empty.places.traps.capacity(), 3);
        assert_eq!(empty.places.holes.capacity(), 6);
        assert_eq!(empty.places.bandtails.capacity(), 9);

        let populated = Cube::new_random_from_density(2.0, 1.0, 1.0, 1.5, 2, 3, true)
            .expect("small random cube should be generated");
        assert_eq!(populated.places.traps.len(), 3);
        assert_eq!(populated.places.holes.len(), 6);
        assert_eq!(populated.places.bandtails.len(), 9);

        let error = Cube::new_empty_from_density(
            1.0,
            1.0,
            1.0,
            2.0,
            usize::MAX,
            0,
            false,
        )
        .unwrap_err();
        assert_eq!(error, "hole count overflowed usize");
    }
}
