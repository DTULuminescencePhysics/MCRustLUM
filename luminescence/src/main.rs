use common::crystal::{Cube, Coord, Boundary};
use std::rc::Rc;

fn main() {
    // Example 1: Padded boundary (default) with random features
    let cube_padded = Cube::new_random(10.0, 10.0, 10.0, 2, 2, 1, false);
    println!("=== Padded Boundary Cube ===");
    println!("Cube dimensions: {}x{}x{}", cube_padded.boundary.x, cube_padded.boundary.y, cube_padded.boundary.z);
    println!("Boundary: {:?}", cube_padded.boundary.kind);
    println!("Traps: {}, Holes: {}, Bandtails: {}\n", 
        cube_padded.places.traps.len(), 
        cube_padded.places.holes.len(),
        cube_padded.places.bandtails.len());

    // Example 2: Periodic boundary with random features
    let cube_periodic = Cube::new_random(10.0, 10.0, 10.0, 2, 2, 1, true);
    println!("=== Periodic Boundary Cube ===");
    println!("Cube dimensions: {}x{}x{}", cube_periodic.boundary.x, cube_periodic.boundary.y, cube_periodic.boundary.z);
    println!("Boundary: {:?}", cube_periodic.boundary.kind);
    println!("Traps: {}, Holes: {}, Bandtails: {}\n", 
        cube_periodic.places.traps.len(), 
        cube_periodic.places.holes.len(),
        cube_periodic.places.bandtails.len());

    // Example 3: Distance calculation with different boundaries
    let p1 = Coord::new(1.0, 5.0, 5.0, cube_periodic.boundary.clone());
    let p2 = Coord::new(9.0, 5.0, 5.0, cube_periodic.boundary.clone());

    println!("=== Distance Comparison ===");
    println!("Point 1: ({}, {}, {})", p1.x, p1.y, p1.z);
    println!("Point 2: ({}, {}, {})", p2.x, p2.y, p2.z);
    
    let padded_distance = cube_padded.distance(&p1, &p2);
    let periodic_distance = cube_periodic.distance(&p1, &p2);

    println!("Distance (Padded):  {:.2}", padded_distance);
    println!("Distance (Periodic): {:.2}", periodic_distance);
    println!("\nWith periodic boundaries, the distance wraps around the box!");

    // Example 4: Using integers for coordinates
    println!("\n=== Integer Coordinates ===");
    let int_boundary = Rc::new(Boundary::new(10i32, 10i32, 10i32, false));
    let int_p1 = Coord::new(1i32, 5i32, 5i32, int_boundary.clone());
    let int_p2 = Coord::new(4i32, 5i32, 5i32, int_boundary);
    let int_distance = int_p1.distance(&int_p2);
    println!("Integer coordinates: ({}, {}, {}) to ({}, {}, {})", 
        int_p1.x, int_p1.y, int_p1.z, int_p2.x, int_p2.y, int_p2.z);
    println!("Distance: {:.2}", int_distance);

    // Example 5: Mixed numeric types (x: i32, y: f32, z: i64)
    println!("\n=== Mixed Numeric Types ===");
    let mixed_boundary = Rc::new(Boundary::new(10i32, 10.0f32, 10i64, false));
    let mixed_p1 = Coord::new(1i32, 5.0f32, 5i64, mixed_boundary.clone());
    let mixed_p2 = Coord::new(4i32, 5.0f32, 5i64, mixed_boundary);
    let mixed_distance = mixed_p1.distance(&mixed_p2);
    println!("Mixed types (i32, f32, i64): ({}, {}, {}) to ({}, {}, {})", 
        mixed_p1.x, mixed_p1.y, mixed_p1.z, mixed_p2.x, mixed_p2.y, mixed_p2.z);
    println!("Distance: {:.2}", mixed_distance);
}
