// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Equations governing the numerical models.
//!
//! This crate contains rate equations. File parsing and Monte Carlo orchestration live in the sibling
//! `io` and `mc` crates.
//! This crate contains crystal geometry, time/temperature profiles, 
//! and the evaluation of the rate equations. File parsing 
//! and Monte Carlo orchestration live in the sibling
//! `io` and `mc` .


/// Crystal geometry and randomly generated electron-site positions.
pub mod crystal;
/// Piecewise-linear time and temperature profiles.
pub mod time_temperature;
/// Spatial sites and the physical parameters assigned to electron traps.
pub mod trap_hole_band_tail;
/// Runtime selection and composition of rate equations.
pub mod rate_equation_selection;
/// Scalar and container-aware physical rate equations.
pub mod rate_equations;
/// Typed parameter groups consumed by rate-equation selections.
pub mod rate_equation_inputs;
/// Charge-transfer events, rate candidates, and sampled lifetimes.
pub mod charge_transfer;


#[cfg(test)]
mod tests {
    use super::*;
    use common::random;
    #[test]
    /// Verify that cube contains random point.
    fn cube_contains_random_point() {
        let cube = crystal::Cube::new(10.0, 20.0, 30.0, vec![0], 0, 0, true).unwrap();
        let mut rng = random::get_std_rng_for_rep(0);
        let point = trap_hole_band_tail::Coord::random_in(10.0,20.0,30.0, &mut rng).unwrap();
        assert!(cube.contains(&point));
    }

    #[test]
    /// Verify that filled random creates features.
    fn filled_random_creates_features() {
        let cube = crystal::Cube::new(5.0, 5.0, 5.0, vec![2], 3, 1, true).unwrap();
        let mut rng = random::get_std_rng_for_rep(0);
        let places = trap_hole_band_tail::ElectronPlaces::random_from_cube(&cube, &mut rng).unwrap();
        assert_eq!(places.traps().len(), 2);
        assert_eq!(places.holes().len(), 3);
        assert_eq!(places.bandtails().unwrap().len(), 1);

        assert!(places
                .traps()
                .iter()
                .all(|trap| cube.contains(&trap)));
        assert!(places
                .holes()
                .iter()
                .all(|hole| cube.contains(&hole)));
        assert!(places
                .bandtails()
                .unwrap()
                .iter()
                .all(|bandtail| cube.contains(&bandtail)));
    }

    #[test]
    /// Verify that mixed types coordinates.
    fn mixed_types_coordinates() {
        // x: i32, y: f32, z: i64
        let boundary = crystal::Boundary::new(10i32, 10.0f32, 10i64, false).unwrap();
        let p1 = trap_hole_band_tail::Coord::new(1i32, 5.0f32, 5i64, ).unwrap();
        let p2 = trap_hole_band_tail::Coord::new(4i32, 5.0f32, 5i64, ).unwrap();

        let distance = boundary.distance(&p1, &p2);
        assert_eq!(distance, 3.0);
    }
}
