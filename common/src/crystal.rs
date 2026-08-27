// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only
//! Spatial representation of the simulation crystal.
//!
//! A [`crate::crystal::Cube`] combines its boundary condition with the
//! number of electron traps, hole recombination sites, and bandtail states.

use crate::numeric::{Float, Numeric};
use crate::trap_hole_band_tail::Coord;

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



/// A simulation volume containing the boundary conditions and number of electron-sites.
///
#[derive(Debug, Clone)]
pub struct Cube {
    /// Trap, hole, and bandtail coordinates in this realization.
    pub trap_total: usize,
    pub hole_total: usize,
    pub bandtail_total: usize,
    /// Dimensions and boundary condition used by distance calculations.
    pub boundary: Boundary,
}

impl Cube {
    /// Create a validated cube.
    pub fn new<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        trap_total: usize,
        hole_total: usize,
        bandtail_total: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let boundary = Boundary::new(x, y, z, periodic)?;
        Ok(Self { trap_total, hole_total, bandtail_total, boundary })
    }
    /// Create an empty cube whose capacities are derived from trap density.
    ///
    /// `h_no` and `b_no` are counts per trap, not absolute counts.
    pub fn new_from_density<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        density: Float,
        h_no: usize,
        b_no: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let boundary = Boundary::new(x, y, z, periodic)?;
        let trap_total = boundary.density_to_number(density)?;
        let hole_total = h_no
            .checked_mul(trap_total)
            .ok_or_else(|| "hole count overflowed usize".to_string())?;
        let bandtail_total = b_no
            .checked_mul(trap_total)
            .ok_or_else(|| "bandtail count overflowed usize".to_string())?;

        Ok(Self { trap_total, hole_total, bandtail_total, boundary })
    }
   
    /// Return whether `point` lies within this cube.
    pub fn contains(&self, point: &Coord) -> bool {
        self.boundary.contains(point)
    }

    /// Calculate the distance between two points using the cube's boundary condition.
    pub fn distance(&self, p1: &Coord, p2: &Coord) -> Float {
        self.boundary.distance(p1, p2)
    }

   
   

}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn coord_with_integers() {
        let boundary = Boundary::new(10i32, 10i32, 10i32, false).unwrap();
        let p1 = Coord::new(0i32, 0i32, 0i32).unwrap();
        let p2 = Coord::new(3i32, 4i32, 0i32).unwrap();
        boundary.distance(&p1, &p2);
        assert_eq!(p1.distance(&p2), 5.0);
    }

    #[test]
    fn generation_rejects_invalid_boundaries_and_density() {
        assert!(Boundary::new(0.0, 10.0, 10.0, false).is_err());
        assert!(Boundary::new(Float::NAN, 10.0, 10.0, false).is_err());
        assert!(Cube::new_from_density(
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

    
}
