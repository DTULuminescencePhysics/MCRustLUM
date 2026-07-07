pub mod numeric;
pub mod crystal;
pub mod time_temperature;
pub mod constants;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cube_contains_random_point() {
        let cube = crystal::Cube::new_random(10.0, 20.0, 30.0, 2, 2, 0, true);
        let point = cube.random_point();
        assert!(cube.contains(&point));
    }

    #[test]
    fn filled_random_creates_features() {
        let cube = crystal::Cube::new_random(5.0, 5.0, 5.0, 2, 3, 1, true);
        assert_eq!(cube.places.traps.len(), 2);
        assert_eq!(cube.places.holes.len(), 3);
        assert_eq!(cube.places.bandtails.len(), 1);
        assert!(cube.places.traps.iter().all(|trap| cube.contains(&trap.position)));
        assert!(cube.places.holes.iter().all(|hole| cube.contains(&hole.position)));
        assert!(cube.places.bandtails.iter().all(|bandtail| cube.contains(&bandtail.position)));
    }

    #[test]
    fn mixed_types_coordinates() {
        // x: i32, y: f32, z: i64
        let boundary = crystal::Boundary::new(10i32, 10.0f32, 10i64, false);
        let p1 = crystal::Coord::new(1i32, 5.0f32, 5i64, );
        let p2 = crystal::Coord::new(4i32, 5.0f32, 5i64, );

        let distance = boundary.distance(&p1, &p2);
        assert_eq!(distance, 3.0);
    }
}
