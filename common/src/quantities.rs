// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

/// Holds structs that control how values are input 
/// and converted between different units
    
use crate::units::{Unit, UnitValue};
use crate::random::{
    BoundedExponentialRandom, GaussianRandom,
    generate_bounded_exponential_random, generate_gaussian_random,
};
use serde::{Deserialize, Deserializer}; 
use rand::Rng;


/// Struct that holds input values 
#[derive(Debug, Clone)]
pub struct Quantity<T = f64> {
    /// Quantity value
    pub value: T,
    /// Unit of the value
    pub unit: Unit,
    /// Flag to decide if the value should be 
    /// generated from a distribution 
    pub distribute: bool,
    /// Optional distribution function to 
    /// generate said values. Cleared during deserialization when `distribute` is false.
    pub function: Option<Distribution>,
}

impl<'de, T: Deserialize<'de>> Deserialize<'de> for Quantity<T> {
    /// Read a quantity and retain its distribution only when explicitly enabled.
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        /// Input representation before applying the distribution switch.
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct QuantityInput<T> {
            value: T,
            unit: Unit,
            #[serde(default)]
            distribute: bool,
            function: Option<Distribution>,
        }

        let input = QuantityInput::<T>::deserialize(deserializer)?;
        Ok(Self {
            value: input.value,
            unit: input.unit,
            distribute: input.distribute,
            function: if input.distribute { input.function } else { None },
        })
    }
}

impl<T: UnitValue> Quantity<T> {

    /// Function to change the unit of a quantity 
    pub fn change_unit(&mut self, new_unit: Unit)-> Result<(), String> {
        use crate::units::Dimension::*;
        let source = self.unit.dimension;
        let target = new_unit.dimension;
        if !self.unit.dimension.same_kind(&new_unit.dimension) {
            return Err(format!("Attempting to convert between {} and {} ",self.unit, new_unit));
        }
        if self.unit == new_unit {
            return Ok(());
        }
        let (scale, offset) = match (source, target) {
            (Temperature{kind}, Temperature{kind:kind_new}) => {
                let offset = kind.convert_unit(&kind_new, 0.0).ok_or("err")?;
                (1.0, offset)
            },

            (Energy{kind}, Energy{kind: kind_new}) => { 
                let scale = kind.convert_unit(kind_new,self.unit.get_scaler() / new_unit.get_scaler()).ok_or("err")?;
                
                (scale,0.0)
            },
            (Frequency{ kind }, Frequency { kind: kind_new } ) => {
                (kind_new.get_time_unit_metric()/kind.get_time_unit_metric(),0.0) 
            },
            (DoseRate{ kind }, DoseRate { kind: kind_new } ) => {
                (kind_new.get_time_unit_metric()/kind.get_time_unit_metric(),0.0) 
            },
            (Time{ kind }, Time { kind: kind_new } ) => {
                (kind.get_time_unit_metric()/kind_new.get_time_unit_metric(),0.0) 
            },
            _ if source == target => {
                let current = self.unit.get_scaler();
                let new = new_unit.get_scaler();

                let scale = match source {
                    InverseLength | InverseArea | Density => new / current,
                    _ => current / new,
                };

                (scale, 0.0)
            }

            _ => {
                return Err(format!(
                    "cannot convert {source:?} to {target:?}"
                ));
            }
        };

        if !scale.is_finite() || scale <= 0.0 {
            println!("{:?}, {:?}", self.unit.dimension, new_unit.dimension);
            return Err(format!("invalid unit conversion scale: {scale}"));
        }

        // Validate before modifying the quantity.
        if self.distribute && self.function.is_none() {
            return Err(
                "distribution selected but no distribution present".into(),
            );
        }
        self.value.apply_conversion(scale, offset);
        if self.distribute {
            match self.function.as_mut() {
                Some(Distribution::Gaussian { sigma }) => {
                    // Standard deviation has no temperature offset.
                    *sigma *= scale;
                }

                Some(Distribution::Exponential { scale: spread, min, max }) => {
                    *min = *min * scale + offset;
                    *max = *max * scale + offset;

                    // The exponential scale has the same units as its bounds.
                    *spread *= scale;
                }

                None => unreachable!("distribution was validated above"),
            }
        }

        self.unit = new_unit;
        Ok(())

    }

}


/// Distribution can either be Gaussian
/// or exponential (further distributions could be added) 
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "kind", deny_unknown_fields)]
pub enum Distribution<T= f64> {
    /// Gaussian distribution which would be centred 
    /// round the [`Quantity`] value
    #[serde(rename = "gaussian")]
    Gaussian { 
        /// Sigma is the width of the gaussian 
        sigma: T 
    },
    
    /// Exponential distribution
    #[serde(rename = "exp", alias = "exponential")]
    Exponential { 
        /// Scale of the bounded exponential distribution.
        #[serde(alias = "rate")]
        scale: T,
        /// minimum value the exponential can take
        min:  T, 
        /// maximum value the exponential can take
        max:  T 
    },
}

impl<T> Distribution<T>
where
    T: Copy + BoundedExponentialRandom + GaussianRandom,
{
    /// Generate a sample using this distribution and the supplied center value.
    pub fn generate_value(&self, value: &T, rng: &mut impl Rng) -> T {
        match self {
            Distribution::Exponential { scale, min, max } => {
                generate_bounded_exponential_random(*min, *max, *scale, rng)
            },
            Distribution::Gaussian { sigma } => {
                generate_gaussian_random(*value, *sigma, rng)
            }
        }


    }

}



#[cfg(test)]
mod tests {
    use super::*;

   #[test]
    /// Verify that changing a temperature quantity converts every array value and its unit.
    fn temperature_quantity_converts_array_values_and_unit() {
        let kelvin: Unit = "kelvin".parse().unwrap();
        let mut temperatures: Quantity<Vec<f64>> = Quantity {
            value: vec![0.0, 800.0],
            unit: "celsius".parse().unwrap(),
            distribute: false,
            function: None,
        };

        temperatures.change_unit(kelvin.clone()).unwrap();
        assert_eq!(temperatures.value, vec![273.15, 1073.15]);
        assert_eq!(temperatures.unit, kelvin);

        // Reapplying normalization must not add the offset a second time.
        temperatures.change_unit(kelvin).unwrap();
        assert_eq!(temperatures.value, vec![273.15, 1073.15]);

        let celsius: Unit = "celsius".parse().unwrap();
        temperatures.change_unit(celsius.clone()).unwrap();
        assert_eq!(temperatures.unit, celsius);
        for (actual, expected) in temperatures.value.iter().zip([0.0, 800.0]) {
            assert!((actual - expected).abs() < 1.0e-12);
        }
    }

    #[test]
    /// Verify that Celsius conversion preserves a scalar's Gaussian spread.
    fn temperature_quantity_preserves_gaussian_spread() {

        let mut temperature = Quantity {
            value: 20.0,
            unit: "celsius".parse().unwrap(),
            distribute: true,
            function: Some(Distribution::Gaussian { sigma: 5.0 }),
        };

        temperature.change_unit("kelvin".parse().unwrap()).unwrap();
        assert_eq!(temperature.value, 293.15);
        assert!(matches!(temperature.function, Some(Distribution::Gaussian { sigma: 5.0 })));
    }

}