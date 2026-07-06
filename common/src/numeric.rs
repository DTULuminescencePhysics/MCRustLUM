use rand::Rng;
/// Trait for numeric types that can be used as coordinates.
/// Supports conversion to f32 for distance calculations and comparisons.
pub trait Numeric: Copy + Clone + PartialOrd + std::fmt::Debug {
    fn to_f32(self) -> f32;
    fn zero() -> Self;
    /// Generate a random value in the range [0, max]
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> f32;
}

impl Numeric for f32 {
    fn to_f32(self) -> f32 { self }
    fn zero() -> Self { 0.0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> f32 {
       rng.gen_range(0.0..=max)
        
    }
}
impl Numeric for f64 {
    fn to_f32(self) -> f32 { self as f32 }
    fn zero() -> Self { 0.0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> f32 {
        let out = rng.gen_range(0.0..=max).to_f32();
        out
    }
}

impl Numeric for i32 {
    fn to_f32(self) -> f32 { self as f32 }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> f32 {
        let out = rng.gen_range(0..=max).to_f32();
        out
    }
}

impl Numeric for i64 {
    fn to_f32(self) -> f32 { self as f32 }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> f32 {
        let out = rng.gen_range(0..=max).to_f32();
        out
    }
}