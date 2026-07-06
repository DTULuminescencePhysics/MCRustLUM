/// This module holds the simulation crystal 

use crate::numeric::Numeric;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryCondition {
    Periodic,
    Padded,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boundary {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub kind: BoundaryCondition,
}

impl Boundary {
    /// Create a Periodic or Padded boundary condition.
    pub fn new<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z, periodic: bool) -> Self {
        let x = x.to_f32();
        let y =y.to_f32();
        let z = z.to_f32();
        
        if periodic {
            Self { x, y, z, kind: BoundaryCondition::Periodic }
        } else {
            Self { x, y, z, kind: BoundaryCondition::Padded }
        }
    }

    /// Padded boundary distance between p1 and p2
    pub fn padded_distance(p1: &Coord, p2: &Coord) -> f32 {
        ((p1.x- p2.x).powi(2) + (p1.y - p2.y).powi(2) + (p1.z - p2.z).powi(2)).sqrt()
    }

    /// Periodic boundary distance between p1 and p2
    pub fn periodic_distance(&self, p1: &Coord, p2: &Coord) -> f32 {
        
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
    pub fn distance(&self, p1: &Coord, p2: &Coord) -> f32 {
        match self.kind {
            BoundaryCondition::Padded => Boundary::padded_distance(p1, p2),
            BoundaryCondition::Periodic => self.periodic_distance(p1, p2),
        }
    }

    pub fn contains(&self, point: &Coord) -> bool {
        point.x >= 0.0 && point.x <= self.x
            && point.y >= 0.0 && point.y <= self.y
            && point.z >= 0.0 && point.z <= self.z
    }
}

/// The coordinate struct holds the x, y, z coordinates and can either be set directly 
/// or randomly generated between some 0 and x,y,z limits
#[derive(Debug, Clone)]
pub struct Coord {
    pub x: f32,
    pub y: f32,
    pub z: f32,
   
}

impl Coord {
    /// Create a Coord 
    pub fn new<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z) -> Self {
        let x = x.to_f32();
        let y =y.to_f32();
        let z = z.to_f32();
        Self { x, y, z }
    }

    /// Randomly generates x, y, z coordinates of a Coord between [0,X), [0,y), [0,z)
    pub fn random_in<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z ) -> Self {
        let mut rng = rand::thread_rng();
        Coord {
            x: X::random_in(x, &mut rng),
            y: Y::random_in(y, &mut rng),
            z: Z::random_in(z, &mut rng),

        }
    }

    /// Calculate distance to another coordinate using this coordinate's boundary.
    pub fn distance(&self, other: &Coord) -> f32 {
        ((self.x- other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)).sqrt()
    }
}

/// Describes an electron trapping site 
#[derive(Debug, Clone)]
pub struct Trap {
    pub position: Coord,
}

impl Trap {
    pub fn new(position: Coord) -> Self {
        Self { position }
    }

    pub fn random<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z,) -> Self {
        let coord = Coord::random_in(x, y, z,);
        Self::new(coord)
    }
}
/// Describes a hole trap capable of accepting an electron 
#[derive(Debug, Clone)]
pub struct Hole {
    pub position: Coord,
}

impl Hole {
    pub fn new(position: Coord) -> Self {
        Self { position }
    }

    pub fn random<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z,) -> Self {
        let coord = Coord::random_in(x, y, z, );
        Self::new(coord)
    }
}

/// Describes a band tail state that are shallow electron traps
#[derive(Debug, Clone)]
pub struct Bandtail {
    pub position: Coord,
}

impl Bandtail {
    pub fn new(position: Coord) -> Self {
        Self { position }
    }

    pub fn random<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z,) -> Self {
        let coord = Coord::random_in(x, y, z, );
        Self::new(coord)
    }
}

/// Struct to store the traps, holes and bandtail states in a simulation crystal
#[derive(Debug, Clone)]
pub struct ElectronPlaces {
    pub traps: Vec<Trap>,
    pub holes: Vec<Hole>,
    pub bandtails: Vec<Bandtail>,
}

impl ElectronPlaces{
    pub fn new(
        traps: Vec<Trap>,
        holes: Vec<Hole>,
        bandtails: Vec<Bandtail>,
    ) -> Self {
        Self { traps, holes, bandtails }
    }

    /// Create a new ElectronPlaces with reserved capacity for later assignment.
    pub fn with_capacity(traps: usize, holes: usize, bandtails: usize) -> Self {
        Self {
            traps: Vec::with_capacity(traps),
            holes: Vec::with_capacity(holes),
            bandtails: Vec::with_capacity(bandtails),
        }
    }

    /// Replace the trap vector later.
    pub fn set_traps(&mut self, traps: Vec<Trap>) {
        self.traps = traps;
    }

    /// Replace the hole vector later.
    pub fn set_holes(&mut self, holes: Vec<Hole>) {
        self.holes = holes;
    }

    /// Replace the bandtail vector later.
    pub fn set_bandtails(&mut self, bandtails: Vec<Bandtail>) {
        self.bandtails = bandtails;
    }

    /// Push a pre-built trap into the collection.
    pub fn push_trap(&mut self, trap: Trap) {
        self.traps.push(trap);
    }

    /// Create and append a trap at the given coordinate.
    pub fn push_trap_at<X: Numeric, Y: Numeric, Z: Numeric>(&mut self, x: X, y: Y, z: Z,) {
        self.traps.push(Trap::new(Coord::new(x, y, z, )));
    }

    /// Push a pre-built hole into the collection.
    pub fn push_hole(&mut self, hole: Hole) {
        self.holes.push(hole);
    }

    /// Create and append a hole at the given coordinate.
    pub fn push_hole_at<X: Numeric, Y: Numeric, Z: Numeric>(&mut self, x: X, y: Y, z: Z, ) {
        self.holes.push(Hole::new(Coord::new(x, y, z, )));
    }

    /// Push a pre-built bandtail into the collection.
    pub fn push_bandtail(&mut self, bandtail: Bandtail) {
        self.bandtails.push(bandtail);
    }

    /// Create and append a bandtail at the given coordinate.
    pub fn push_bandtail_at<X: Numeric, Y: Numeric, Z: Numeric>(&mut self, x: X, y: Y, z: Z, ) {
        self.bandtails.push(Bandtail::new(Coord::new(x, y, z, )));
    }

    /// Calculate distance between two traps.
    pub fn trap_trap_distance(&self, p1: usize, p2: usize) -> f32 {
        self.traps[p1].position.distance(&self.traps[p2].position)
    }

    /// Create a new ElectronPlaces with random coordinates.
    pub fn random_new<X: Numeric, Y: Numeric, Z: Numeric>(
        t_no: usize,
        h_no: usize,
        b_no: usize,
        x: X,
        y: Y,
        z: Z,
    
    ) -> Self {
        let traps = (0..t_no)
            .map(|_| Trap::random(x, y, z,))
            .collect();

        let holes = (0..h_no)
            .map(|_| Hole::random(x, y, z, ))
            .collect();

        let bandtails = (0..b_no)
            .map(|_| Bandtail::random(x, y, z,))
            .collect();

        ElectronPlaces::new(traps, holes, bandtails)
    }
}

#[derive(Debug, Clone)]
pub struct Cube {
    pub places: ElectronPlaces,
    pub boundary: Boundary,
}

impl Cube {
    pub fn new<X: Numeric, Y: Numeric, Z: Numeric>(x: X, y: Y, z: Z, periodic: bool) -> Self {
        let boundary = Boundary::new(x, y, z, periodic);
        let places = ElectronPlaces::with_capacity(0, 0, 0);
        Self { places, boundary }
    }
    pub fn new_random<X: Numeric, Y: Numeric, Z: Numeric>(
        x: X,
        y: Y,
        z: Z,
        t_no: usize,
        h_no: usize,
        b_no: usize,
        periodic: bool,
    ) -> Self {
        let boundary = Boundary::new(x, y, z, periodic);
        let places = ElectronPlaces::random_new(t_no, h_no, b_no, x, y, z,);
        Self { places, boundary }
    }
    pub fn contains(&self, point: &Coord) -> bool {
        self.boundary.contains(point)
    }

    /// Calculate distance between two points with the cube's boundary conditions
    pub fn distance(&self, p1: &Coord, p2: &Coord) -> f32 {
        self.boundary.distance(p1, p2)
    }

    pub fn trap_trap_distance(&self, p1: usize, p2: usize) -> f32 {
        self.boundary.distance(&self.places.traps[p1].position,&self.places.traps[p2].position)
    }

    pub fn random_point(&self) -> Coord {
        Coord::random_in(self.boundary.x, self.boundary.y, self.boundary.z, )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_euclidean() {
        let p1 = Coord::new(0.0, 0.0, 0.0,);
        let p2 = Coord::new(3.0, 4.0, 0.0,);
        assert_eq!(p1.distance(&p2), 5.0);
    }

    #[test]
    fn distance_padded_boundary() {
        let p1 = Coord::new(1.0, 1.0, 1.0,);
        let p2 = Coord::new(4.0, 1.0, 1.0,);
        let boundary = Boundary::new(10.0, 10.0, 10.0, false);
        let distance = boundary.distance(&p1, &p2);
        assert_eq!(distance, 3.0);
    }

    #[test]
    fn distance_periodic_boundary_wraps() {
        let p1 = Coord::new(0.5, 5.0, 5.0,);
        let p2 = Coord::new(8.5, 5.0, 5.0,);
        let boundary = Boundary::new(10.0, 10.0, 10.0, true);
        // Without wrapping: distance = 9.0
        // With wrapping: distance should be 2.0 (shorter path wraps around)
       let distance = boundary.distance(&p1, &p2);
        assert_eq!(distance, 2.0);
    }

    #[test]
    fn cube_distance_with_boundary() {
        let mut cube = Cube::new_random(10.0, 10.0, 10.0, 0, 0, 0, true);
        cube.places.push_trap_at(1.0, 5.0, 5.0,);
        cube.places.push_trap_at(9.0, 5.0, 5.0,);
        // Periodic: shortest distance is 2.0
        assert_eq!(cube.trap_trap_distance(0, 1), 2.0);
    }

    #[test]
    fn coord_with_integers() {
        let boundary = Boundary::new(10i32, 10i32, 10i32, false); 
        let p1 = Coord::new(0i32, 0i32, 0i32, );
        let p2 = Coord::new(3i32, 4i32, 0i32,);
        boundary.distance(&p1, &p2);
        assert_eq!(p1.distance(&p2), 5.0);
    }
}