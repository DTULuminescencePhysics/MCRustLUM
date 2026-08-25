//! Numeric types and element-wise operations shared by the simulation.
//!
//! Rate equations are written once and can then operate on a scalar, a
//! [`Vec`], or an [`ndarray::ArrayBase`]. [`crate::numeric::ElementWise`]
//! supplies the binary arithmetic used by those equations, while
//! [`crate::numeric::ElementWiseUnary`] supplies
//! operations such as `exp` and `powf`.
//!
//! Two floating-point precisions are intentionally available:
//! [`crate::numeric::Float`] is the normal simulation precision and
//! [`crate::numeric::TimeFloat`] is the higher precision used for rates and
//! time calculations. [`crate::numeric::PrecisionInput`] converts either
//! individual values or complete containers between them.

use std::ops::{Add, Div, Mul, Sub};
use rand::Rng;
use ndarray::{Array, Array1, ArrayBase, ArrayD, Data, Dimension, Zip};

/// Default floating-point precision used for model parameters and state.
pub type Float = f64;
/// Higher precision used for rates, lifetimes, and time/temperature profiles.
pub type TimeFloat = f64;

/// Common behaviour required from scalar values accepted by numeric helpers.
///
/// This trait allows coordinates and model inputs to be supplied as common
/// integer or floating-point types. Conversion is explicit so calculations
/// can consistently target either [`Float`] or [`TimeFloat`].
pub trait Numeric: Copy + Clone + PartialOrd + std::fmt::Debug {
    /// Convert this value to the normal simulation precision.
    fn to_float(self) -> Float;
    /// Convert this value to the precision used for time and rates.
    fn to_time_float(self) -> TimeFloat;
    /// Return the additive identity for this numeric type.
    fn zero() -> Self;
    /// Generate a value in the inclusive range `0..=max`.
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

/// Describes a precision to which numeric inputs can be converted.
///
/// Marker types implement this trait so [`PrecisionInput`] can select its
/// output type without duplicating conversion logic for scalars and
/// containers.
pub trait PrecisionTarget {
    /// Scalar value produced by this precision target.
    type Value: Copy + Mul<Output = Self::Value>;
    /// Convert one supported scalar into the target precision.
    fn from_numeric<T: Numeric>(value: T) -> Self::Value;

}

/// Marker selecting the program-wide [`Float`] precision.
pub struct ProgramPrecision;
/// Marker selecting the higher-precision [`TimeFloat`] representation.
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

/// Converts a scalar or container to a selected [`PrecisionTarget`].
///
/// Implementations preserve the input shape: a scalar remains a scalar, a
/// vector remains a vector, and an ndarray retains its dimensions. Both
/// methods consume the input because conversion creates a new value.
pub trait PrecisionInput<P>
where
    P: PrecisionTarget,
{
    /// Converted scalar or shape-preserving container type.
    type Output;
    /// Convert every value to the selected precision.
    fn map_to_precision(self) -> Self::Output;
    /// Convert every value and multiply it by a target-precision scalar.
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

/// Floating-point operations needed by generic rate equations.
///
/// Unlike [`Numeric`], this trait is restricted to real floating-point types
/// because exponential and fractional-power operations are required.
pub trait RealNumber:
    Numeric
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
{
    /// Calculate `e` raised to this value.
    fn exp(self) -> Self;
    /// Raise this value to a floating-point power.
    fn powf(self, power: Self) -> Self;
    /// Construct this type from an `f64`, accepting a precision loss for f32.
    fn from_f64(value: f64) -> Self;
    /// Convert any supported numeric scalar to `f64`.
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

/// Shape-aware binary arithmetic for scalars, vectors, and ndarrays.
///
/// Scalar/container combinations broadcast the scalar across the container.
/// Two vectors must have equal lengths and two ndarrays must have identical
/// shapes. Vector/ndarray combinations use ndarray broadcasting rules and
/// return a dynamically dimensioned array. Operations return `None` when
/// their input shapes cannot be combined; arithmetic conditions such as
/// division by zero follow the underlying floating-point behaviour.
pub trait ElementWise<Rhs = Self> {
    /// Result type after applying the operation and any broadcasting.
    type Output;

    /// Add corresponding elements.
    fn element_add(&self, rhs: &Rhs) -> Option<Self::Output>;
    /// Subtract corresponding elements, preserving left-to-right order.
    fn element_sub(&self, rhs: &Rhs) -> Option<Self::Output>;
    /// Multiply corresponding elements.
    fn element_mul(&self, rhs: &Rhs) -> Option<Self::Output>;
    /// Divide corresponding elements, preserving left-to-right order.
    fn element_div(&self, rhs: &Rhs) -> Option<Self::Output>;
}   

/// Unary mathematical operations applied independently to every element.
///
/// Unary operations cannot have shape mismatches, so they return their output
/// directly rather than wrapping it in [`Option`].
pub trait ElementWiseUnary {
    /// Shape-preserving result of the unary operation.
    type Output;

    /// Apply the exponential function to every value.
    fn element_exp(&self) -> Self::Output;
    /// Raise every value to the same numeric power.
    fn element_powf<P>(&self, power: P) -> Self::Output
    where
        P: Numeric;
}

// Scalar operations cannot fail because there is no container shape to check.
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

// Mixed f32/f64 arithmetic promotes the result to f64 so the more precise
// operand is not narrowed before the operation.
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

impl ElementWise<Vec<f64>> for f32 {
    type Output = Vec<f64>;

    fn element_add(&self, rhs: &Vec<f64>) -> Option<Self::Output> {
        let left = *self as f64;
        Some(rhs.iter().map(|value| left + *value).collect())
    }

    fn element_sub(&self, rhs: &Vec<f64>) -> Option<Self::Output> {
        let left = *self as f64;
        Some(rhs.iter().map(|value| left - *value).collect())
    }

    fn element_mul(&self, rhs: &Vec<f64>) -> Option<Self::Output> {
        let left = *self as f64;
        Some(rhs.iter().map(|value| left * *value).collect())
    }

    fn element_div(&self, rhs: &Vec<f64>) -> Option<Self::Output> {
        let left = *self as f64;
        Some(rhs.iter().map(|value| left / *value).collect())
    }
}

impl ElementWise<f32> for Vec<f64> {
    type Output = Vec<f64>;

    fn element_add(&self, rhs: &f32) -> Option<Self::Output> {
        let right = *rhs as f64;
        Some(self.iter().map(|value| *value + right).collect())
    }

    fn element_sub(&self, rhs: &f32) -> Option<Self::Output> {
        let right = *rhs as f64;
        Some(self.iter().map(|value| *value - right).collect())
    }

    fn element_mul(&self, rhs: &f32) -> Option<Self::Output> {
        let right = *rhs as f64;
        Some(self.iter().map(|value| *value * right).collect())
    }

    fn element_div(&self, rhs: &f32) -> Option<Self::Output> {
        let right = *rhs as f64;
        Some(self.iter().map(|value| *value / right).collect())
    }
}

// Paired vectors require equal lengths. Using `Option` prevents `zip` from
// silently discarding values from the longer input.
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

// A scalar is broadcast over every vector element in either operand order.
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

// Array/array operations deliberately require exact shapes. More permissive
// broadcasting is reserved for the explicit vector/array implementations
// below, making accidental array shape changes less likely.
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

// Scalar/array operations preserve the ndarray's static dimension type.
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

// Combining Vec and ndarray uses ndarray broadcasting. Because the final
// dimensionality depends on the array supplied at runtime, the output is
// ArrayD rather than an array with a statically known dimension.
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

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn precision_conversion_preserves_scalar_vector_and_array_shapes() {
        let scalar = <i32 as PrecisionInput<ProgramPrecision>>::map_to_precision(3);
        let vector = <Vec<i16> as PrecisionInput<TimePrecision>>::multiply_to_precision(
            vec![1, 2, 3],
            2.5,
        );
        let matrix = <ndarray::Array2<i32> as PrecisionInput<ProgramPrecision>>::map_to_precision(
            array![[1, 2], [3, 4]],
        );

        assert_eq!(scalar, 3.0);
        assert_eq!(vector, vec![2.5, 5.0, 7.5]);
        assert_eq!(matrix, array![[1.0, 2.0], [3.0, 4.0]]);
    }

    #[test]
    fn scalar_vector_operations_preserve_operand_order() {
        let vector = vec![2.0_f64, 4.0];

        assert_eq!(10.0_f32.element_sub(&vector), Some(vec![8.0, 6.0]));
        assert_eq!(10.0_f32.element_div(&vector), Some(vec![5.0, 2.5]));
        assert_eq!(vector.element_sub(&2.0_f32), Some(vec![0.0, 2.0]));
        assert_eq!(vector.element_div(&2.0_f32), Some(vec![1.0, 2.0]));
    }

    #[test]
    fn paired_containers_reject_incompatible_shapes() {
        let short = vec![1.0_f64, 2.0];
        let long = vec![1.0_f64, 2.0, 3.0];
        let matrix = array![[1.0_f64, 2.0], [3.0, 4.0]];
        let row = array![[1.0_f64, 2.0]];

        assert!(short.element_add(&long).is_none());
        assert!(matrix.element_mul(&row).is_none());
    }

    #[test]
    fn vector_array_operations_broadcast_the_trailing_dimension() {
        let vector = vec![10.0_f64, 20.0, 30.0];
        let matrix = array![[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let sum = vector.element_add(&matrix).expect("vector should broadcast");
        let difference = matrix.element_sub(&vector).expect("vector should broadcast");

        assert_eq!(sum, array![[11.0, 22.0, 33.0], [14.0, 25.0, 36.0]].into_dyn());
        assert_eq!(
            difference,
            array![[-9.0, -18.0, -27.0], [-6.0, -15.0, -24.0]].into_dyn(),
        );
        assert!(vec![1.0_f64, 2.0].element_add(&matrix).is_none());
    }

    #[test]
    fn unary_operations_apply_to_every_container_element() {
        let vector = vec![1.0_f64, 2.0, 3.0];
        let squared = vector.element_powf(2_i32);
        let exponentials = array![0.0_f64, 1.0].element_exp();

        assert_eq!(squared, vec![1.0, 4.0, 9.0]);
        assert!((exponentials[0] - 1.0).abs() < 1.0e-12);
        assert!((exponentials[1] - std::f64::consts::E).abs() < 1.0e-12);
    }
}
