use std::ops::{Add, Div, Mul, Sub};
use rand::Rng;
use ndarray::{Array, Array1, ArrayBase, ArrayD, Data, Dimension, Zip};
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
        value.to_time_float()
    }
}

pub trait PrecisionInput<P>
where
    P: PrecisionTarget,
{
    type Output;
    fn map_to_precision(self) -> Self::Output;
    fn multiply_to_precision(self, multiplier: P::Value) -> Self::Output;
}

impl<T, P> PrecisionInput<P> for T
where
    T: Numeric,
    P: PrecisionTarget,
{
    type Output = P::Value;

    fn map_to_precision(self) -> Self::Output {
        P::from_numeric(self)
    }

    fn multiply_to_precision(self, multiplier: P::Value) -> Self::Output {
        P::from_numeric(self) * multiplier
    }
}

impl<T, P> PrecisionInput<P> for Vec<T>
where
    T: Numeric,
    P: PrecisionTarget,
{
    type Output = Vec<P::Value>;

    fn map_to_precision(self) -> Self::Output {
        self.into_iter()
            .map(P::from_numeric)
            .collect()
    }

    fn multiply_to_precision(self, multiplier: P::Value) -> Self::Output {
        self.into_iter()
            .map(|value| P::from_numeric(value) * multiplier)
            .collect()
    }
}

impl<T, P, S, D> PrecisionInput<P> for ArrayBase<S, D>
where
    T: Numeric,
    P: PrecisionTarget,
    S: Data<Elem = T>,
    D: Dimension,
{
    type Output = Array<P::Value, D>;

    fn map_to_precision(self) -> Self::Output {
        self.mapv(P::from_numeric)
    }

    fn multiply_to_precision(self, multiplier: P::Value) -> Self::Output {
        self.mapv(|value| P::from_numeric(value) * multiplier)
    }
}

pub trait RealNumber:
    Numeric
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
{
    fn exp(self) -> Self;
    fn powf(self, power: Self) -> Self;
    fn from_f64(value: f64) -> Self;
    fn to_f64<T>(value:T)-> f64 
    where
        T: Numeric ;
}

impl RealNumber for f32 {
    fn exp(self) -> Self {
        f32::exp(self)
    }

    fn powf(self, power: Self) -> Self {
        f32::powf(self, power)
    }

    fn from_f64(value: f64) -> Self {
        value as f32
    }
    fn to_f64<T>(value: T) -> f64 
    where 
        T:Numeric,
    {
        value.to_float() as f64
    }
}

impl RealNumber for f64 {
    fn exp(self) -> Self {
        f64::exp(self)
    }

    fn powf(self, power: Self) -> Self {
        f64::powf(self, power)
    }

    fn from_f64(value: f64) -> Self {
        value
    }

    fn to_f64<T>(value: T) -> Self 
    where 
        T:Numeric,
    {
        value.to_float() as Self
    }

}

pub trait ElementWise<Rhs = Self> {
    type Output;

    fn element_add(&self, rhs: &Rhs) -> Option<Self::Output>;
    fn element_sub(&self, rhs: &Rhs) -> Option<Self::Output>;
    fn element_mul(&self, rhs: &Rhs) -> Option<Self::Output>;
    fn element_div(&self, rhs: &Rhs) -> Option<Self::Output>;
}   

pub trait ElementWiseUnary {
    type Output;

    fn element_exp(&self) -> Self::Output;
    fn element_powf<P>(&self, power: P) -> Self::Output
    where
        P: Numeric;
}

macro_rules! impl_element_wise_scalar {
    ($type:ty) => {
        impl ElementWise<$type> for $type {
            type Output = $type;

            fn element_add(&self, rhs: &$type) -> Option<Self::Output> {
                Some(*self + *rhs)
            }

            fn element_sub(&self, rhs: &$type) -> Option<Self::Output> {
                Some(*self - *rhs)
            }

            fn element_mul(&self, rhs: &$type) -> Option<Self::Output> {
                Some(*self * *rhs)
            }

            fn element_div(&self, rhs: &$type) -> Option<Self::Output> {
                Some(*self / *rhs)
            }
        }

        impl ElementWiseUnary for $type {
            type Output = $type;

            fn element_exp(&self) -> Self::Output {
                <$type>::exp(*self)
            }

            fn element_powf<P>(&self, power: P) -> Self::Output 
            where 
                P: Numeric,
            {   
                let power = <$type as RealNumber>::to_f64(power);
                <$type>::powf(*self, power as $type)
            }
        }
    };
}

impl_element_wise_scalar!(f32);
impl_element_wise_scalar!(f64);

impl ElementWise<f32> for f64 {
    type Output = f64;

    fn element_add(&self, rhs: &f32) -> Option<Self::Output> {
        Some(*self + *rhs as f64)
    }

    fn element_sub(&self, rhs: &f32) -> Option<Self::Output> {
        Some(*self - *rhs as f64)
    }

    fn element_mul(&self, rhs: &f32) -> Option<Self::Output> {
        Some(*self * *rhs as f64)
    }

    fn element_div(&self, rhs: &f32) -> Option<Self::Output> {
        Some(*self / *rhs as f64)
    }
}

impl ElementWise<f64> for f32 {
    type Output = f64;

    fn element_add(&self, rhs: &f64) -> Option<Self::Output> {
        Some(*self as f64 + *rhs)
    }

    fn element_sub(&self, rhs: &f64) -> Option<Self::Output> {
        Some(*self as f64 - *rhs)
    }

    fn element_mul(&self, rhs: &f64) -> Option<Self::Output> {
        Some(*self as f64 * *rhs)
    }

    fn element_div(&self, rhs: &f64) -> Option<Self::Output> {
        Some(*self as f64 / *rhs)
    }
}

impl<T> ElementWise<Vec<T>> for Vec<T>
where
    T: RealNumber,
{
    type Output = Vec<T>;

    fn element_add(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        (self.len() == rhs.len()).then(|| {
            self.iter()
                .zip(rhs.iter())
                .map(|(left, right)| *left + *right)
                .collect()
        })
    }

    fn element_sub(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        (self.len() == rhs.len()).then(|| {
            self.iter()
                .zip(rhs.iter())
                .map(|(left, right)| *left - *right)
                .collect()
        })
    }

    fn element_mul(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        (self.len() == rhs.len()).then(|| {
            self.iter()
                .zip(rhs.iter())
                .map(|(left, right)| *left * *right)
                .collect()
        })
    }

    fn element_div(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        (self.len() == rhs.len()).then(|| {
            self.iter()
                .zip(rhs.iter())
                .map(|(left, right)| *left / *right)
                .collect()
        })
    }
}

impl<T> ElementWise<T> for Vec<T>
where
    T: RealNumber,
{
    type Output = Vec<T>;

    fn element_add(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.iter().map(|value| *value + *rhs).collect())
    }

    fn element_sub(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.iter().map(|value| *value - *rhs).collect())
    }

    fn element_mul(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.iter().map(|value| *value * *rhs).collect())
    }

    fn element_div(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.iter().map(|value| *value / *rhs).collect())
    }
}

impl<T> ElementWise<Vec<T>> for T
where
    T: RealNumber,
{
    type Output = Vec<T>;

    fn element_add(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        Some(rhs.iter().map(|value| *self + *value).collect())
    }

    fn element_sub(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        Some(rhs.iter().map(|value| *self - *value).collect())
    }

    fn element_mul(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        Some(rhs.iter().map(|value| *self * *value).collect())
    }

    fn element_div(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        Some(rhs.iter().map(|value| *self / *value).collect())
    }
}

impl<T> ElementWiseUnary for Vec<T>
where
    T: RealNumber,
{
    type Output = Vec<T>;

    fn element_exp(&self) -> Self::Output {
        self.iter().map(|value| value.exp()).collect()
    }

    fn element_powf<P>(&self, power:P) -> Self::Output 
    where 
        P:Numeric,
    {   
        let power = T::from_f64(T::to_f64(power));
        self.iter().map(|value| value.powf(power)).collect()
    }
}

impl<T, S, D, S2> ElementWise<ArrayBase<S2, D>> for ArrayBase<S, D>
where
    T: RealNumber,
    S: Data<Elem = T>,
    S2: Data<Elem = T>,
    D: Dimension,
{
    type Output = Array<T, D>;

    fn element_add(&self, rhs: &ArrayBase<S2, D>) -> Option<Self::Output> {
        (self.shape() == rhs.shape()).then(|| Zip::from(self).and(rhs).map_collect(|left, right| *left + *right))
    }

    fn element_sub(&self, rhs: &ArrayBase<S2, D>) -> Option<Self::Output> {
        (self.shape() == rhs.shape()).then(|| Zip::from(self).and(rhs).map_collect(|left, right| *left - *right))
    }

    fn element_mul(&self, rhs: &ArrayBase<S2, D>) -> Option<Self::Output> {
        (self.shape() == rhs.shape()).then(|| Zip::from(self).and(rhs).map_collect(|left, right| *left * *right))
    }

    fn element_div(&self, rhs: &ArrayBase<S2, D>) -> Option<Self::Output> {
        (self.shape() == rhs.shape()).then(|| Zip::from(self).and(rhs).map_collect(|left, right| *left / *right))
    }
}

impl<T, S, D> ElementWise<T> for ArrayBase<S, D>
where
    T: RealNumber,
    S: Data<Elem = T>,
    D: Dimension,
{
    type Output = Array<T, D>;

    fn element_add(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.mapv(|value| value + *rhs))
    }

    fn element_sub(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.mapv(|value| value - *rhs))
    }

    fn element_mul(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.mapv(|value| value * *rhs))
    }

    fn element_div(&self, rhs: &T) -> Option<Self::Output> {
        Some(self.mapv(|value| value / *rhs))
    }
}

impl<T, S, D> ElementWise<ArrayBase<S, D>> for T
where
    T: RealNumber,
    S: Data<Elem = T>,
    D: Dimension,
{
    type Output = Array<T, D>;

    fn element_add(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        Some(rhs.mapv(|value| *self + value))
    }

    fn element_sub(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        Some(rhs.mapv(|value| *self - value))
    }

    fn element_mul(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        Some(rhs.mapv(|value| *self * value))
    }

    fn element_div(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        Some(rhs.mapv(|value| *self / value))
    }
}

impl<T, S, D> ElementWiseUnary for ArrayBase<S, D>
where
    T: RealNumber,
    S: Data<Elem = T>,
    D: Dimension,
{
    type Output = Array<T, D>;

    fn element_exp(&self) -> Self::Output {
        self.mapv(|value| value.exp())
    }

    fn element_powf<P>(&self, power: P) -> Self::Output 
    where 
        P:Numeric,
    {
        let power = T::from_f64(T::to_f64(power));
        self.mapv(|value| value.powf(power))
    }
}

impl<T, S, D> ElementWise<ArrayBase<S, D>> for Vec<T>
where
    T: RealNumber,
    S: Data<Elem = T>,
    D: Dimension,
{
    type Output = ArrayD<T>;

    fn element_add(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        let left = Array1::from(self.clone()).into_dyn();
        let right = rhs.view().into_dyn();
        let left = left.broadcast(right.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left + *right))
    }

    fn element_sub(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        let left = Array1::from(self.clone()).into_dyn();
        let right = rhs.view().into_dyn();
        let left = left.broadcast(right.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left - *right))
    }

    fn element_mul(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        let left = Array1::from(self.clone()).into_dyn();
        let right = rhs.view().into_dyn();
        let left = left.broadcast(right.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left * *right))
    }

    fn element_div(&self, rhs: &ArrayBase<S, D>) -> Option<Self::Output> {
        let left = Array1::from(self.clone()).into_dyn();
        let right = rhs.view().into_dyn();
        let left = left.broadcast(right.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left / *right))
    }
}

impl<T, S, D> ElementWise<Vec<T>> for ArrayBase<S, D>
where
    T: RealNumber,
    S: Data<Elem = T>,
    D: Dimension,
{
    type Output = ArrayD<T>;

    fn element_add(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        let left = self.view().into_dyn();
        let right = Array1::from(rhs.clone()).into_dyn();
        let right = right.broadcast(left.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left + *right))
    }

    fn element_sub(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        let left = self.view().into_dyn();
        let right = Array1::from(rhs.clone()).into_dyn();
        let right = right.broadcast(left.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left - *right))
    }

    fn element_mul(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        let left = self.view().into_dyn();
        let right = Array1::from(rhs.clone()).into_dyn();
        let right = right.broadcast(left.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left * *right))
    }

    fn element_div(&self, rhs: &Vec<T>) -> Option<Self::Output> {
        let left = self.view().into_dyn();
        let right = Array1::from(rhs.clone()).into_dyn();
        let right = right.broadcast(left.raw_dim())?;
        Some(Zip::from(left).and(right).map_collect(|left, right| *left / *right))
    }
}
