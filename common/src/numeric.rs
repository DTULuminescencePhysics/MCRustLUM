use std::ops::Mul;
use rand::Rng;
/// Trait for numeric types that can be used as coordinates.
/// Set the default precision used by most of the program.
pub type Float = f32;

/// Set the precision used by time/temperature profiles.
pub type TimeFloat = f64;

// pub type Float = f64;
/// Supports conversion to the program-wide `Float` type for distance calculations and comparisons.
pub trait Numeric: Copy + Clone + PartialOrd + std::fmt::Debug {
    fn to_float(self) -> Float;
    fn to_time_float(self) -> TimeFloat;
    fn zero() -> Self;
    /// Generate a random value in the range [0, max]
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float;
}

impl Numeric for f32 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0.0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
       let out = rng.gen_range(0.0..=max).to_float();
       out
    }
}
impl Numeric for f64 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0.0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0.0..=max).to_float();
        out
    }
}
impl Numeric for i8 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for i16 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for i32 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for i64 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for u8 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for u16 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for u32 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

impl Numeric for u64 {
    fn to_float(self) -> Float { self as Float }
    fn to_time_float(self) -> TimeFloat { self as TimeFloat }
    fn zero() -> Self { 0 }
    fn random_in(max: Self, rng: &mut rand::rngs::ThreadRng) -> Float {
        let out = rng.gen_range(0..=max).to_float();
        out
    }
}

pub trait PrecisionTarget {
    type Value: Copy + Mul<Output = Self::Value>;
    fn from_numeric<T: Numeric>(value: T) -> Self::Value;

}

pub struct ProgramPrecision;
pub struct TimePrecision;

impl PrecisionTarget for ProgramPrecision {
    type Value = Float;
    fn from_numeric<T: Numeric>(value: T) -> Self::Value {
        value.to_float()
    }
}

impl PrecisionTarget for TimePrecision {
    type Value = TimeFloat;
    fn from_numeric<T: Numeric>(value: T) -> Self::Value {
        value.to_float() as TimeFloat
    }
}

pub trait PrecisionInput<P>
where
    P: PrecisionTarget,
{
    type Output;
    fn map_to_precision(self, multiplier: P::Value) -> Self::Output;
}

impl<T, P> PrecisionInput<P> for T
where
    T: Numeric,
    P: PrecisionTarget,
{
    type Output = P::Value;

    fn map_to_precision(self, multiplier: P::Value) -> Self::Output {
        P::from_numeric(self) * multiplier
    }
}

impl<T, P> PrecisionInput<P> for Vec<T>
where
    T: Numeric,
    P: PrecisionTarget,
{
    type Output = Vec<P::Value>;

    fn map_to_precision(self, multiplier: P::Value) -> Self::Output {
        self.into_iter()
            .map(|value| P::from_numeric(value) * multiplier)
            .collect()
    }
}
