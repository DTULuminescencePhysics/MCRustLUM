// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! SI conversions and physical constants used by the transport equations.
/// Internal scalar, square, and cubic SI prefix factors.
///
/// These constants are currently reserved for future unit-aware geometry.
pub mod metric{
    use std::{fmt:: Debug, str::FromStr};    

    /// Factor represented by the SI prefix giga.
    const GIGA: f64 = 1_000_000_000.0;
    /// Factor represented by the SI prefix mega.
    const MEGA: f64 = 1_000_000.0;
    /// Factor represented by the SI prefix kilo.
    const KILO: f64 = 1_000.0;
    /// Factor represented by the SI prefix hecto.
    const HECTO: f64 = 100.0;
    /// Factor represented by the SI prefix deca.
    const DECA: f64 = 10.0;
    /// Factor represented by the SI prefix deci.
    const DECI: f64 = 0.1;
    /// Factor represented by the SI prefix centi.
    const CENTI: f64 = 0.01;
    /// Factor represented by the SI prefix milli.
    const MILLI: f64 = 0.001;
    /// Factor represented by the SI prefix micro.
    const MICRO: f64 = 0.000_001;
    /// Factor represented by the SI prefix nano.
    const NANO: f64 = 0.000_000_001;
    /// Length factor for one ångström in metres.
    const ANG: f64 = 0.0_000_000_001;
    /// Factor represented by the SI prefix pico.
    const PICO: f64 = 0.000_000_000_001;
    /// Factor represented by the SI prefix femto.
    const FEMTO: f64 = 0.000_000_000_000_001;

    /// Squared giga factor for area conversions.
    const GIGA_SQ: f64 = GIGA*GIGA;
    /// Squared mega factor for area conversions.
    const MEGA_SQ: f64 = MEGA*MEGA;
    /// Squared kilo factor for area conversions.
    const KILO_SQ: f64 = KILO*KILO;
    /// Squared hecto factor for area conversions.
    const HECTO_SQ: f64 = HECTO*HECTO;
    /// Squared deca factor for area conversions.
    const DECA_SQ: f64 = DECA*DECA;
    /// Squared deci factor for area conversions.
    const DECI_SQ: f64 = DECI*DECI;
    /// Squared centi factor for area conversions.
    const CENTI_SQ: f64 = CENTI*CENTI;
    /// Squared milli factor for area conversions.
    const MILLI_SQ: f64 = MILLI*MILLI;
    /// Squared micro factor for area conversions.
    const MICRO_SQ: f64 = MICRO*MICRO;
    /// Squared nano factor for area conversions.
    const NANO_SQ: f64 = NANO*NANO;
    /// Squared ångström factor for area conversions.
    const ANG_SQ: f64 = ANG*ANG;
    /// Squared pico factor for area conversions.
    const PICO_SQ: f64 = PICO*PICO;
    /// Squared femto factor for area conversions.
    const FEMTO_SQ: f64 = FEMTO*FEMTO;

    /// Cubed mega factor for volume conversions.
    const MEGA_CUBE: f64 = MEGA_SQ*MEGA;
    /// Cubed kilo factor for volume conversions.
    const KILO_CUBE: f64 = KILO_SQ*KILO;
    /// Cubed hecto factor for volume conversions.
    const HECTO_CUBE: f64 = HECTO_SQ*HECTO;
    /// Cubed deca factor for volume conversions.
    const DECA_CUBE: f64 = DECA_SQ*DECA;
    /// Cubed deci factor for volume conversions.
    const DECI_CUBE: f64 = DECI_SQ*DECI;
    /// Cubed centi factor for volume conversions.
    const CENTI_CUBE: f64 = CENTI_SQ*CENTI;
    /// Cubed milli factor for volume conversions.
    const MILLI_CUBE: f64 = MILLI_SQ*MILLI;
    /// Cubed micro factor for volume conversions.
    const MICRO_CUBE: f64 = MICRO_SQ*MICRO;
    /// Cubed nano factor for volume conversions.
    const NANO_CUBE: f64 = NANO_SQ*NANO;
    /// Cubed ångström factor for volume conversions.
    const ANG_CUBE: f64 = ANG_SQ*ANG;
    /// Cubed pico factor for volume conversions.
    const PICO_CUBE: f64 = PICO_SQ*PICO;
    /// Cubed femto factor for volume conversions.
    const FEMTO_CUBE: f64 = FEMTO_SQ*FEMTO;

    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum metric_abbrev {
        GIGA,
        MEGA,
        KILO,
        HECTO,
        DECA,
        DECI,
        CENTI,
        STANDARD,
        MILLI,
        MICRO,
        NANO,
        ANG,
        PICO,
        FEMTO,
    }
     impl FromStr for metric_abbrev {
        type Err = String;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s.to_lowercase().as_str() {
                "giga"   => Ok(Self::GIGA),
                "mega"   => Ok(Self::MEGA),
                "kilo"   => Ok(Self::KILO),
                "hecto"  => Ok(Self::HECTO),
                "deca"  => Ok(Self::DECA),
                "deci" => Ok(Self::DECI),
                "milli" => Ok(Self::MILLI),
                "micro" => Ok(Self::MICRO),
                "nano" => Ok(Self::NANO),
                "ang" => Ok(Self::ANG),
                "pico" => Ok(Self::PICO),
                "femto" => Ok(Self::FEMTO),
                _ => Ok(Self::STANDARD),
            }
        }
    }


    impl metric_abbrev {

        pub fn get_name(&self)-> &'static str {
            match self {
                metric_abbrev::STANDARD => "",
                metric_abbrev::GIGA => "G",
                metric_abbrev::MEGA => "M",
                metric_abbrev::KILO => "k",
                metric_abbrev::HECTO => "h",
                metric_abbrev::DECA => "da",
                metric_abbrev::DECI => "d",
                metric_abbrev::CENTI => "c",
                metric_abbrev::MILLI => "m",
                metric_abbrev::MICRO => "μ",
                metric_abbrev::NANO => "n",
                metric_abbrev::ANG => "Å",
                metric_abbrev::PICO => "p",
                metric_abbrev::FEMTO => "f",
            }
        }
        pub fn get_metric(&self)-> f64 {
            match self {
                metric_abbrev::STANDARD => 1.0,
                metric_abbrev::GIGA => GIGA,
                metric_abbrev::MEGA => MEGA,
                metric_abbrev::KILO => KILO,
                metric_abbrev::HECTO => HECTO,
                metric_abbrev::DECA => DECA,
                metric_abbrev::DECI => DECI,
                metric_abbrev::CENTI => CENTI,
                metric_abbrev::MILLI => MILLI,
                metric_abbrev::MICRO => MICRO,
                metric_abbrev::NANO => NANO,
                metric_abbrev::ANG => ANG,
                metric_abbrev::PICO => PICO,
                metric_abbrev::FEMTO => FEMTO,
            }
        }
        pub fn get_square_metric(&self)-> f64 {
            match self {
                metric_abbrev::STANDARD => 1.0,
                metric_abbrev::GIGA => GIGA_SQ,
                metric_abbrev::MEGA => MEGA_SQ,
                metric_abbrev::KILO => KILO_SQ,
                metric_abbrev::HECTO => HECTO_SQ,
                metric_abbrev::DECA => DECA_SQ,
                metric_abbrev::DECI => DECI_SQ,
                metric_abbrev::CENTI => CENTI_SQ,
                metric_abbrev::MILLI => MILLI_SQ,
                metric_abbrev::MICRO => MICRO_SQ,
                metric_abbrev::NANO => NANO_SQ,
                metric_abbrev::ANG => ANG_SQ,
                metric_abbrev::PICO => PICO_SQ,
                metric_abbrev::FEMTO => FEMTO_SQ,
            }
        }
        pub fn get_cube_metric(&self)-> f64 {
            match self {
                metric_abbrev::STANDARD => 1.0,
                metric_abbrev::GIGA => 1.0,
                metric_abbrev::MEGA => MEGA_CUBE,
                metric_abbrev::KILO => KILO_CUBE,
                metric_abbrev::HECTO => HECTO_CUBE,
                metric_abbrev::DECA => DECA_CUBE,
                metric_abbrev::DECI => DECI_CUBE,
                metric_abbrev::CENTI => CENTI_CUBE,
                metric_abbrev::MILLI => MILLI_CUBE,
                metric_abbrev::MICRO => MICRO_CUBE,
                metric_abbrev::NANO => NANO_CUBE,
                metric_abbrev::ANG => ANG_CUBE,
                metric_abbrev::PICO => PICO_CUBE,
                metric_abbrev::FEMTO => FEMTO_CUBE,
            }
        }

    }

}

/// Fundamental constants expressed in SI units, plus Boltzmann's constant in eV/K.
pub mod physical_constants{
    use crate::numeric::{Float};
    /// Boltzmann constant in J/K (kg·m²/s²·K)
    pub const BOLTZMANN: Float = 1.380649e-23;
    /// Boltzmann constant in eV/K
    pub const BOLTZMANN_EV: Float = 8.617333262e-5;
    /// Planck constant in J·s (kg·m²/s)
    pub const PLANCK: Float = 6.62607015e-34;
    /// Reduced Planck constant in J·s (kg·m²/s)
    pub const PLANCK_BAR: Float = 1.054571817e-34;
    /// Speed of light in m/s
    pub const SPEED_LIGHT: Float = 299792458.0;
    /// Elementary charge in C
    pub const ELEMENTARY_CHARGE: Float = 1.602176634e-19;
    /// Electron mass in kg
    pub const ELECTRON_MASS: Float = 9.10938356e-31;
    /// Avogadro's number in 1/mol
    pub const AVOGARDRO: Float = 6.02214076e23;
    /// Gas constant in J/(mol·K)
    pub const GAS_CONSTANT: Float = 8.314462618;
    /// Vacuum permittivity in F/m
    pub const VACUUM_PERM: Float = 8.854187817e-12;
    /// Conversion factor from eV to J
    pub const EV_TO_J: Float = 1.602176634e-19;
}

/// Time-unit conversion at the precision used by the simulation clock.
pub mod time{
    use std::{fmt::{self, Debug}, str::FromStr};    
    use crate::numeric::{ElementWise, Float, PrecisionInput, TimeFloat, TimePrecision,};
    /// Seconds in one second.
    const SECOND: f64 = 1.0;
    /// Seconds in one minute.
    const MINUTE: f64 = 60.0;
    /// Seconds in one hour.
    const HOUR: f64 = 60.0*MINUTE;
    /// Seconds in one 24-hour day.
    const DAY: f64 = 24.0*HOUR;
    /// Seconds in one Julian year of 365.2425 days.
    const YEAR: f64 = 31_556_952.0;
    /// Seconds in one thousand Julian years.
    const K_ANNUM: f64 = 31_556_952_000.0;
    /// Seconds in one million Julian years.
    const MA_ANNUM: f64 = 31_556_952_000_000.0;

     /// Unit used for time values supplied to the simulation.
    ///
    /// Geological ages use annum-based names: `KAnnum` is one thousand years
    /// and `MaAnnum` is one million years.
    #[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize, Eq)]
    #[serde(rename_all = "snake_case")]
    pub enum TimeUnit {
        /// Seconds.
        Second,
        /// Minutes.
        Minute,
        /// Hours.
        Hour,
        /// Days of 24 hours.
        Day,
        /// Julian years of 365.2425 days.
        Year,
        /// Thousands of Julian years (`ka`).
        KAnnum,
        /// Millions of Julian years (`Ma`).
        MaAnnum
    }
    impl FromStr for TimeUnit {
        type Err = String;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s.to_lowercase().as_str() {
                "second"   => Ok(Self::Second),
                "minute"   => Ok(Self::Minute),
                "hour"     => Ok(Self::Hour),
                "day"      => Ok(Self::Day),
                "year"     => Ok(Self::Year),
                "ka annum" => Ok(Self::KAnnum),
                "ma annum" => Ok(Self::MaAnnum),
                "kaannum" => Ok(Self::KAnnum),
                "maannum" => Ok(Self::MaAnnum),
                _ => Err(format!("Invalid TimeUnit: {s}")),
            }
        }
    }
    impl fmt::Display for TimeUnit {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            let name = match self {
                TimeUnit::Second => "second",
                TimeUnit::Minute => "minute",
                TimeUnit::Hour => "hour",
                TimeUnit::Day => "day",
                TimeUnit::Year => "year",
                TimeUnit::KAnnum => "Ka",
                TimeUnit::MaAnnum => "Ma",
            };

            formatter.write_str(name)
        }
    }
    impl TimeUnit{
        /// Checks if the time unit is Ka or Ma and so ordering is reversed
        pub fn is_ka_or_ma(&self) -> bool {
            match self {
                TimeUnit::KAnnum  => true,
                TimeUnit::MaAnnum  => true,
                _ => false,
            }
        }
        /// Return the number of seconds in one `unit`.
        pub fn get_time_unit_metric(&self,) -> f64 {
            match self {
                TimeUnit::Second => SECOND,
                TimeUnit::Minute => MINUTE,
                TimeUnit::Hour => HOUR,
                TimeUnit::Day => DAY,
                TimeUnit::Year => YEAR,
                TimeUnit::KAnnum => K_ANNUM,
                TimeUnit::MaAnnum => MA_ANNUM,
            }
        }
        pub fn get_name(&self) -> &'static str {
            match self {
               TimeUnit::Second => "second",
                TimeUnit::Minute => "minute",
                TimeUnit::Hour => "hour",
                TimeUnit::Day => "day",
                TimeUnit::Year => "year",
                TimeUnit::KAnnum => "Ka",
                TimeUnit::MaAnnum => "Ma",
            }
        }
        pub fn convert_unit<T>(&self, unit: &TimeUnit, value:T)-> Option<T::Output>
        where
            T: PrecisionInput<TimePrecision>,
            T::Output: ElementWise<Float, Output = T::Output>,
        {
            let value = value.map_to_precision();
            if self == unit {
                Some(value)
            } else{
                let mult = unit.get_time_unit_metric()/self.get_time_unit_metric();
                value.element_mul(&mult)
            }
        }  

    }
     /// Return the number of seconds in one `unit`.
    pub fn get_time_unit_metric(name: &str) -> f64 {
        if name == "second" {
            SECOND
        }else if name == "minute"{
            MINUTE
        }else if name == "hour"{
            HOUR
        }else if name == "day"{
            DAY
        }else if name == "year"{
            YEAR
        }else if name == "Ka"{
            K_ANNUM
        }else if name == "Ma"{
            MA_ANNUM
        } else {
            SECOND
        }
       
    }
    /// Return the number of seconds in one `unit`.
    pub fn unit_multiplier(unit: TimeUnit) -> TimeFloat {
        get_time_unit_metric(&unit.to_string()) as TimeFloat
    }

    /// Convert a scalar or supported container to seconds at time precision.
    ///
    /// The same function accepts scalar and vector values:
    ///
    /// ```
    /// use common::constants::time::{convert_to_seconds, TimeUnit};
    ///
    /// let seconds = convert_to_seconds(TimeUnit::Minute, vec![1.0, 2.0]);
    /// assert_eq!(seconds, Some(vec![60.0, 120.0]));
    /// ```
    pub fn convert_to_seconds<T>(unit: TimeUnit, value: T) -> Option<T::Output>
    where
        T: PrecisionInput<TimePrecision>,
    {
        Some(value.multiply_to_precision(unit_multiplier(unit)))
    }
}

/// Temperature-unit conversion at the program's configured precision.
pub mod temperature {
    use std::{fmt, str::FromStr};        
    use crate::numeric::{ElementWise, Float, PrecisionInput, ProgramPrecision,};
    /// Celsius to kelvin conversion value
    const CELSIUS_TO_KELVIN: Float = 273.15;

    #[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize, Eq)]
    #[serde(rename_all = "snake_case")]
    /// Unit used for temperature values supplied to the simulation.
    pub enum TemperatureUnit {
        /// Degrees Celsius.
        Celsius,
        /// Kelvin.
        Kelvin,
    }
    impl fmt::Display for TemperatureUnit {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            let name = match self {
                TemperatureUnit::Celsius => "°C",
                TemperatureUnit::Kelvin => "K",
            };

            formatter.write_str(name)
        }
    }
    impl FromStr for TemperatureUnit {
        type Err = String;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s.to_lowercase().as_str() {
                "kelvin" | "k" => Ok(Self::Kelvin),
                "celsius" | "°c" => Ok(Self::Celsius),
                _ => Err(format!("Invalid TemperatureUnit: {s}")),
            }
        }
    }
    impl TemperatureUnit {

        pub fn convert_unit<T>(&self, unit: &TemperatureUnit, value:T)-> Option<T::Output>
        where
            T: PrecisionInput<ProgramPrecision>,
            T::Output: ElementWise<Float, Output = T::Output>,
        {

            let value = value.map_to_precision();
            match (self, unit) {     
                (TemperatureUnit::Celsius, TemperatureUnit::Kelvin) => {
                    value.element_add(&CELSIUS_TO_KELVIN)
                },
                (TemperatureUnit::Kelvin, TemperatureUnit::Celsius) => {
                    value.element_sub(&CELSIUS_TO_KELVIN)
                },
                (TemperatureUnit::Celsius, TemperatureUnit::Celsius) => Some(value),
                (TemperatureUnit::Kelvin, TemperatureUnit::Kelvin) =>  Some(value),
            
            
            } 

        }   

        pub fn get_name(&self) -> &'static str {
            match self {
                TemperatureUnit::Celsius => "°C",
                TemperatureUnit::Kelvin => "K",
            }
        }
    }


}

/// Holds structs that control how values are input 
/// and converted between different units
pub mod quantities{
    use crate::numeric::{ElementWise, Float, PrecisionInput, ProgramPrecision,};
    use crate::constants::time;
    use crate::constants::temperature;
    use crate::constants::temperature::TemperatureUnit;
    use crate::constants::time::TimeUnit;
    use crate::random::{
        BoundedExponentialRandom, GaussianRandom,
        generate_bounded_exponential_random, generate_gaussian_random,
    };
    use crate::constants::metric::metric_abbrev;
    use serde::{Deserialize, Deserializer};
    use std::{fmt, str::FromStr};        
    use rand::Rng;

    /// A scalar or container whose values support unit conversion.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum EnergyUnit {
        Joule,
        ElectronVolt,
    }
    impl EnergyUnit {
        pub fn get_name(&self,) -> &'static str {
            match self {
                EnergyUnit::Joule => "J",
                EnergyUnit::ElectronVolt => "eV",
               
            }
        }
        
        pub fn convert_unit<T>(&self, unit: EnergyUnit, value:T)-> Option<T::Output>
        where
            T: PrecisionInput<ProgramPrecision>,
            T::Output: ElementWise<Float, Output = T::Output>,
        {
            let value = value.map_to_precision();
            match (self, unit) {     
                (EnergyUnit::ElectronVolt, EnergyUnit::Joule) => {
                    value.element_mul(&crate::constants::physical_constants::EV_TO_J)
                },
                (EnergyUnit::Joule, EnergyUnit::ElectronVolt) => {
                    value.element_div(&crate::constants::physical_constants::EV_TO_J)
                },
                (EnergyUnit::Joule, EnergyUnit::Joule) => Some(value),
                (EnergyUnit::ElectronVolt, EnergyUnit::ElectronVolt) =>  Some(value),
            } 
        }
    }
    
    pub trait UnitValue {
        fn apply_conversion(&mut self, scale: f64, offset: f64);
    }

    impl UnitValue for f64 {
        fn apply_conversion(&mut self, scale: f64, offset: f64) {
            *self = *self * scale + offset;
        }
    }
    impl UnitValue for f32 {
        fn apply_conversion(&mut self, scale: f64, offset: f64) {
            *self = (f64::from(*self) * scale + offset) as f32;
        }
    }

    impl<T: UnitValue> UnitValue for Vec<T> {
        fn apply_conversion(&mut self, scale: f64, offset: f64) {
            for value in self {
                value.apply_conversion(scale, offset);
            }
        }
    }

    impl<T: UnitValue, const N: usize> UnitValue for [T; N] {
        fn apply_conversion(&mut self, scale: f64, offset: f64) {
            for value in self {
                value.apply_conversion(scale, offset);
            }
        }
    }

    // Also support owned ndarrays and mutable ndarray views.
    impl<T, S, D> UnitValue for ndarray::ArrayBase<S, D>
    where
        T: UnitValue,
        S: ndarray::DataMut<Elem = T>,
        D: ndarray::Dimension,
    {
        fn apply_conversion(&mut self, scale: f64, offset: f64) {
            for value in self.iter_mut() {
                value.apply_conversion(scale, offset);
            }
        }
    }

    /// Explicit definition of 
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Dimension {
        /// No dimension 
        Dimensionless,
        /// A length is some scaled metre
        Length,
        /// Reciprocal of a length
        InverseLength,
        /// Area in some scaled square metre
        Area,
        /// Volume in some scaled cubic metre
        Volume,
        /// Reciprocal of an area
        InverseArea,
        /// Number of something per volume
        Density,
        /// Time in any unit amount
        Time {kind: time::TimeUnit},
        /// Reciprocal of time
        Frequency {kind: time::TimeUnit},
        /// Temperature in celsius or kelvin
        Temperature { kind: temperature::TemperatureUnit},
        /// Energy 
        Energy{ kind: EnergyUnit},
        /// Dose in some scaled number of Gy
        Dose,
        /// Dose per time
        DoseRate {kind: time::TimeUnit},
    }

    /// Holds information about the unit and how it is scaled
    #[derive(Debug, Clone, PartialEq)]
    pub struct Unit {
        /// Holds the metric type milli, mega etc. 
        /// or the time unit 
        pub metric: metric_abbrev,
        /// Holds the type of dimension the unit is
        pub dimension: Dimension, 
    }

    impl FromStr for Unit {
        type Err = String;

        fn from_str(input: &str) -> Result<Self, Self::Err> {

            // Normalize capitalization and whitespace.
            let normalised = input
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ")
                .to_ascii_lowercase();

            match normalised.as_str() {
                "" | "1" => Ok(Self {
                                metric: metric_abbrev::STANDARD,
                                dimension: Dimension::Dimensionless,
                            }),
                "second" => Ok(Self {
                                metric: metric_abbrev::STANDARD,
                                dimension: Dimension::Time {kind: time::TimeUnit::Second},
                                }),
                "minute" => Ok(Self {
                                metric: metric_abbrev::STANDARD,
                                dimension: Dimension::Time {kind: time::TimeUnit::Minute},
                            }),
                "hour" => Ok(Self {
                                metric: metric_abbrev::STANDARD,
                                dimension: Dimension::Time {kind: time::TimeUnit::Hour},
                        }),
                "day" => Ok(Self {
                                    metric: metric_abbrev::STANDARD,
                                    dimension: Dimension::Time {kind: time::TimeUnit::Day},
                         }),
                "year"     => Ok(Self {
                            metric: metric_abbrev::STANDARD,
                            dimension: Dimension::Time {kind: time::TimeUnit::Year},
                        }),
                "ka annum" | "kaannum" | "ka" => Ok(Self {
                            metric: metric_abbrev::STANDARD,
                            dimension: Dimension::Time {kind: time::TimeUnit::KAnnum},
                                }),
                "ma annum" | "maannum" | "ma" => Ok(Self {
                            metric: metric_abbrev::STANDARD,
                            dimension: Dimension::Time {kind: time::TimeUnit::MaAnnum},
                        }),
                "kelvin" | "k" => Ok(Self {
                            metric: metric_abbrev::STANDARD,
                            dimension: Dimension::Temperature {kind: temperature::TemperatureUnit::Kelvin},
                        }),
                "celsius" | "c" => Ok(Self {
                            metric: metric_abbrev::STANDARD,
                            dimension: Dimension::Temperature{kind: temperature::TemperatureUnit::Kelvin},
                        }),
                _ => {
                        let words: Vec<&str> = normalised.split_whitespace().collect();
                    
                        let mut contains_slash = None;
                        let mut key_unit_word = None;
                        let mut metric = None;
                        let mut power = None;
                        let mut time = None;
                        // Reciprocal units must start with "1 /".
                        let (reciprocal, words) = match words.as_slice() {
                            ["1", "/", rest @ ..] => (true, rest),
                            words => (false, words),
                        };
                        for &word in words {
                            match word {
                                "/" if contains_slash.is_none() => {
                                    contains_slash = Some(true)
                                }
                                "metre" if key_unit_word.is_none() => {
                                    key_unit_word = Some("metre")
                                }
                                "ev" if key_unit_word.is_none() => {
                                    key_unit_word = Some("ev")
                                }
                                "joule" if key_unit_word.is_none() => {
                                    key_unit_word = Some("joule")
                                }
                              
                                "gy" if key_unit_word.is_none() => {
                                    key_unit_word = Some("gy")
                                }
                                "second" if key_unit_word.is_none() => {
                                    key_unit_word = Some("second")
                                }
                                "minute" if key_unit_word.is_none() => {
                                    key_unit_word = Some("minute")
                                }
                                "hour" if key_unit_word.is_none() => {
                                    key_unit_word = Some("hour")
                                }
                                "day" if key_unit_word.is_none() => {
                                    key_unit_word = Some("day")
                                }
                                "year" if key_unit_word.is_none() => {
                                    key_unit_word = Some("year")
                                }
                                "kaannum" | "ka_annum" | "ka" if key_unit_word.is_none() => {
                                    key_unit_word = Some("Ka")
                                }
                                "maannum" | "ma_annum" | "ma" if key_unit_word.is_none() => {
                                    key_unit_word = Some("Ma")
                                }

                                "second" if time.is_none() => {
                                    time = Some("second")
                                }
                                "minute" if time.is_none() => {
                                    time = Some("minute")
                                }
                                "hour" if time.is_none() => {
                                    time = Some("hour")
                                }
                                "day" if time.is_none() => {
                                    time = Some("day")
                                }
                                "year" if time.is_none() => {
                                    time = Some("year")
                                }
                                "kaannum" | "ka_annum" | "ka" if time.is_none() => {
                                    time = Some("Ka")
                                }
                                "maannum" | "ma_annum" | "ma" if time.is_none() => {
                                    time = Some("Ma")
                                }

                                "square" if power.is_none() => {
                                    power = Some(2);
                                }

                                "cube" if power.is_none() => {
                                    power = Some(3);
                                }

                                "giga" | "mega" | "kilo" | "hecto"
                                | "deca" | "deci" | "centi" | "milli"
                                | "micro" | "nano" | "ang" | "pico"
                                | "femto" if metric.is_none() =>
                                {
                                    metric = Some(word);
                                }

                                _ => {
                                    return Err(format!(
                                        "unexpected or repeated unit component {word:?}"
                                    ));
                                }
                            }
                        }
                        let key_unit_word = key_unit_word.ok_or("Key word should have been found")?;
                        let power = power.unwrap_or(1);
                        let contains_slash = contains_slash.unwrap_or(false);
                        let time = time.unwrap_or("");
                        let metric = metric_abbrev::from_str(metric.unwrap_or(""))?;

                        let dimension =  if key_unit_word == "metre" {
                            match (reciprocal, power) {
                                (false, 1) => Dimension::Length,
                                (false, 2) => Dimension::Area,
                                (false, 3) => Dimension::Volume,
                                (true, 1) => Dimension::InverseLength,
                                (true, 2) => Dimension::InverseArea,
                                (true, 3) => Dimension::Density,
                                _ => unreachable!(),
                            }
                        } else if key_unit_word == "ev" {
                            Dimension::Energy{kind: EnergyUnit::ElectronVolt}
                        } else if key_unit_word == "joule" {
                            Dimension::Energy{kind: EnergyUnit::Joule}
                        } else if key_unit_word == "gy" {
                            if contains_slash{
                                Dimension::Dose
                            }else {
                                Dimension::DoseRate{ kind: time::TimeUnit::from_str(time)?}
                            }
                        }else if key_unit_word == "celsius" {
                            return Ok(Self {
                                metric: metric_abbrev::STANDARD,
                                dimension: Dimension::Temperature {kind: temperature::TemperatureUnit::Celsius},
                            })
                        } else if key_unit_word == "kelvin" {
                            return Ok(Self {
                                metric: metric_abbrev::STANDARD,
                                dimension: Dimension::Temperature {kind: temperature::TemperatureUnit::Kelvin},
                            })
                        } else {
                            if reciprocal {
                                Dimension::Frequency { kind: time::TimeUnit::from_str(key_unit_word)? }

                            } else {
                                Dimension::Time { kind: time::TimeUnit::from_str(key_unit_word)? }
                            }
                        };
                        Ok(Unit {
                            metric,
                            dimension,
                        })  
                    } 
                }}
    }

    impl Unit {

        fn get_scaler(&self) -> f64 {
            match self.dimension {
                Dimension::Dimensionless => 1.0,
                Dimension::Length => self.metric.get_metric(),
                Dimension::InverseLength => self.metric.get_metric(),
                Dimension::Area => self.metric.get_square_metric(),
                Dimension::InverseArea => self.metric.get_square_metric(),
                Dimension::Volume => self.metric.get_cube_metric(),
                Dimension::Density => self.metric.get_cube_metric(),
                Dimension::Time{kind} => kind.get_time_unit_metric(),
                Dimension::Frequency{kind} => kind.get_time_unit_metric(),
                Dimension::DoseRate{kind} => kind.get_time_unit_metric(),
                Dimension::Energy { .. } => self.metric.get_metric(),
                Dimension::Temperature { .. } => self.metric.get_metric(),
                Dimension::Dose => self.metric.get_metric(),
            }
        }
        pub fn get_time_unit(&self) -> Result<TimeUnit,String> {

            match self.dimension {

                Dimension::Time { kind } 
                | Dimension::Frequency { kind } 
                | Dimension::DoseRate { kind } => {
                    Ok(kind)
                },
                _ => Err(format!("No TimeUnit in {:?}",self.dimension))

            }
        } 
        pub fn get_temperature_unit(&self) -> Result<TemperatureUnit,String>{
             match self.dimension {

                Dimension::Temperature { kind }  => {
                    Ok(kind)
                },
                _ => Err(format!("No TemperatureUnit in {:?}",self.dimension))

            }

        }


    }
    impl<'de> Deserialize<'de> for Unit {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            let input = String::deserialize(deserializer)?;

            input.parse().map_err(serde::de::Error::custom)
        }
    }
    impl fmt::Display for Unit {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {

             let (name, name2,name3)= match self.dimension {
                Dimension::Dimensionless => ("","",""),
                Dimension::Length => {
                    (self.metric.get_name(), "m","")
                },
                Dimension::InverseLength => {
                    (self.metric.get_name(), "m","⁻¹")
                },
                Dimension::Area => {
                    (self.metric.get_name(), "m","²")
                },
                Dimension::InverseArea => {
                    (self.metric.get_name(), "m","⁻²")
                },
                Dimension::Volume => {
                    (self.metric.get_name(), "m","³")
                },
                Dimension::Density => {
                    (self.metric.get_name(), "m","⁻³")
                },
                Dimension::Time{kind} => {
                    (kind.get_name(), "","")
                },
                Dimension::Frequency{kind} => {
                    (kind.get_name(), "⁻¹","")
                },
                Dimension::DoseRate{kind} => {
                    ("Gy",kind.get_name(), "⁻¹")
                },
                Dimension::Energy { kind } => {
                    (self.metric.get_name(),kind.get_name(), "")
                },
                Dimension::Temperature { kind } => {
                    (kind.get_name(),"","")
                },
                Dimension::Dose =>{
                    ("Gy", "", "")
                },
            };

            let name = match self.metric {
                metric_abbrev::ANG => {
                    format!("{name}{name3}")
                },
                _ => {
                    format!("{name}{name2}{name3}")
                }
            };
            
            formatter.write_str(&name)
        }
    }

    


    /// Struct that holds input values 
    #[derive(Debug, Clone, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Quantity<T = f64> {
        /// Quantity value
        pub value: T,
        /// Unit of the value
        pub unit: Unit,
        /// Flag to decide if the value should be 
        /// generated from a distribution 
        #[serde(default)]
        pub distribute: bool,
        /// Optional distribution function to 
        /// generate said values
        pub function: Option<Distribution>,
    }

    impl<T: UnitValue> Quantity<T> {

        /// Function to change the unit of a quantity 
        pub fn change_unit(&mut self, new_unit: Unit)-> Result<(), String> {
            use Dimension::*;
            let source = self.unit.dimension;
            let target = new_unit.dimension;
            if self.unit.dimension != new_unit.dimension {
                return Err(format!("Attempting to convert between {} and {} ",self.unit, new_unit));
            }

            let (scale, offset) = match (source, target) {
                (Temperature{kind}, Temperature{kind:kind_new}) => {
                    let offset = kind.convert_unit(&kind_new, 0.0).ok_or("err")?;
                    (0.0, offset)
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

}


#[cfg(test)]
mod tests {
    use super::temperature::TemperatureUnit;
    use super::time::{ TimeUnit, convert_to_seconds, unit_multiplier};
    use crate::numeric::Float;

    #[test]
    fn converts_celsius_to_kelvin_at_program_precision() {
        let temperature = TemperatureUnit::Celsius.convert_unit(&TemperatureUnit::Kelvin, 20.0_f32).unwrap();

        assert_eq!(temperature, 293.15 as Float);
    }

    #[test]
    fn converts_temperature_vectors_to_kelvin() {
        let temperatures = TemperatureUnit::Celsius.convert_unit(&TemperatureUnit::Kelvin, vec![0_i32, 100_i32]).unwrap();


        assert_eq!(temperatures, vec![273.15 as Float, 373.15 as Float]);
    }

    #[test]
    fn preserves_kelvin_values_while_converting_precision() {
        let temperature = TemperatureUnit::Kelvin.convert_unit(&TemperatureUnit::Kelvin, 300.0_f32).unwrap();


        assert_eq!(temperature, 300.0 as Float);
    }

    #[test]
    fn converts_kelvin_to_celsius_for_scalars_and_vectors() {
        let scalar = TemperatureUnit::Kelvin.convert_unit(&TemperatureUnit::Celsius, 373.15_f64).unwrap();
        let vector = TemperatureUnit::Kelvin.convert_unit(&TemperatureUnit::Celsius, vec![273.15_f64, 293.15, 373.15],).unwrap();

       

        assert!((scalar - 100.0).abs() < 1.0e-12);
        assert_eq!(vector, vec![0.0, 20.0, 100.0]);
        assert_eq!(
            TemperatureUnit::Celsius.convert_unit(&TemperatureUnit::Celsius, 25_i32),
            Some(25.0),
        );
    }

    #[test]
    fn temperature_units_parse_and_display_consistently() {
        assert_eq!("kelvin".parse(), Ok(TemperatureUnit::Kelvin));
        assert_eq!("celsius".parse(), Ok(TemperatureUnit::Celsius));
        assert_eq!(TemperatureUnit::Kelvin.to_string(), "K");
        assert_eq!(TemperatureUnit::Celsius.to_string(), "°C");
    }

    #[test]
    fn time_units_display_with_human_readable_names() {
        let names = [
            (TimeUnit::Second, "second"),
            (TimeUnit::Minute, "minute"),
            (TimeUnit::Hour, "hour"),
            (TimeUnit::Day, "day"),
            (TimeUnit::Year, "year"),
            (TimeUnit::KAnnum, "Ka"),
            (TimeUnit::MaAnnum, "Ma"),
        ];

        for (unit, expected) in names {
            assert_eq!(unit.to_string(), expected);
        }
    }

    #[test]
    fn time_multipliers_and_seconds_cover_short_and_geological_units() {
        assert!(matches!(
            unit_multiplier(TimeUnit::Hour),
            3_600.0,
        ));
        assert!(matches!(
            unit_multiplier(TimeUnit::KAnnum),
            31_556_952_000.0,
        ));
        assert!(!(TimeUnit::Year.is_ka_or_ma()));
        assert!((TimeUnit::KAnnum.is_ka_or_ma()));
        assert!((TimeUnit::MaAnnum.is_ka_or_ma()));

        assert_eq!(
            convert_to_seconds(TimeUnit::Day, vec![0_i32, 1, 2]),
            Some(vec![0.0, 86_400.0, 172_800.0]),
        );
        assert_eq!(
            convert_to_seconds(TimeUnit::MaAnnum, 1.0_f64),
            Some(31_556_952_000_000.0),
        );
    }
}
