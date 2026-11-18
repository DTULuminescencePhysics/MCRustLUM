// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use crate::constants::{time::TimeUnit, temperature::TemperatureUnit, physical_constants::EnergyUnit};
use crate::constants::metric::MetricAbbrev;
use serde::{Deserialize, Deserializer};
use std::{fmt, str::FromStr};        

/// Physical dimension and, where needed, its concrete time, temperature, or energy unit.
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
    Time {
        /// Concrete unit used within this physical dimension.
        kind: TimeUnit,
    },
    /// Reciprocal of time
    Frequency {
        /// Concrete unit used within this physical dimension.
        kind: TimeUnit,
    },
    /// Temperature in celsius or kelvin
    Temperature {
        /// Concrete unit used within this physical dimension.
        kind: TemperatureUnit,
    },
    /// Energy 
    Energy{
        /// Concrete unit used within this physical dimension.
        kind: EnergyUnit,
    },
    /// Dose in some scaled number of Gy
    Dose,
    /// Dose per time
    DoseRate {
        /// Concrete unit used within this physical dimension.
        kind: TimeUnit,
    },
}
impl Dimension {
    /// Return whether dimensions share a variant, ignoring their embedded unit choices.
    pub fn same_kind(&self, other: &Self) -> bool {
        std::mem::discriminant(self) == std::mem::discriminant(other)
    }
}

/// Holds information about the unit and how it is scaled
#[derive(Debug, Clone, PartialEq)]
pub struct Unit {
    /// Holds the metric type milli, mega etc. 
    /// or the time unit 
    pub metric: MetricAbbrev,
    /// Holds the type of dimension the unit is
    pub dimension: Dimension, 
}

impl<'de> Deserialize<'de> for Unit {
    /// Deserialize a unit string through the shared parser and report invalid expressions.
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let input = String::deserialize(deserializer)?;

        input.parse().map_err(serde::de::Error::custom)
    }
}

impl FromStr for Unit {
    type Err = String;

    /// Parse a physical-unit expression after normalizing whitespace and ASCII case.
    ///
    /// Recognizes dimensionless values, time/temperature/energy units, metric
    /// prefixes, and supported inverse/dose-rate forms. Unsupported expressions
    /// return an error.
    fn from_str(input: &str) -> Result<Self, Self::Err> {

        let normalised = input.to_ascii_lowercase();
        
        let (numerator, denominator) = match normalised.split_once('/') {
                None => (normalised.trim(), None),
                Some((left, right)) => {
                    let left = left.trim();
                    let right = right.trim();
                    if right.contains('/') {
                        return Err(format!("Multiple '/' separators in {input:?}"));
                    }
                    if left.is_empty() || right.is_empty() {
                        return Err(format!(
                            "Expected a unit on both sides of '/' in {input:?}"
                        ));
                    }
                    (left, Some(right))
                }
        };
        let Some(denominator) = denominator else {
            return Self::parse_single_unit(numerator);
        };

        let top = Self::parse_single_unit(numerator)?;
        let bottom = Self::parse_single_unit(denominator)?;

        let (metric, dimension) = match (top.dimension, bottom.dimension) {
            (Dimension::Dimensionless, Dimension::Length) => (bottom.metric, Dimension::InverseLength),
            (Dimension::Dimensionless, Dimension::Area) => (bottom.metric, Dimension::InverseArea),
            (Dimension::Dimensionless, Dimension::Volume) => (bottom.metric, Dimension::Density),
            (Dimension::Dimensionless, Dimension::Time { kind }) => (top.metric, Dimension::Frequency { kind }),
            (Dimension::Dose, Dimension::Time { kind }) => (top.metric, Dimension::DoseRate { kind }),
            _ => return Err(format!("Unsupported unit division in {input:?}")),
        };
       
        Ok(Self { metric, dimension })

    }
}
impl Unit {
    /// Helper function used by `[Unit::from_str]` to parse input values either side of 
    /// the slash symbol
    fn parse_single_unit(input: &str) -> Result<Self, String> {
        
        let words: Vec<&str> = input.split_whitespace().collect();
        if words.is_empty() || words.as_slice() == ["1"] {
            return Ok(Self {
                metric: MetricAbbrev::STANDARD,
                dimension: Dimension::Dimensionless,
            });
        }
        let words = match words.as_slice() {
            ["ka", "annum"] => vec!["ka"],
            ["ma", "annum"] => vec!["ma"],
            _ => words,
        };
        
        let mut distance: bool = false;
        let mut radiation: bool = false;
        let mut power: Option<u8> = None;
        let mut time: Option< &str > = None;
        let mut metric: Option< &str > = None;
        let mut energy: Option< &str > = None;
        let mut temperature: Option< &str > = None;
        for word in words {
            match word {
                "1" => {
                    return Err(format!("If a 1 is specified it cannot have any other unit word with it to make it clear \n 
                                        it is a reciprocal unit instead: \n
                                        {input}. \n 
                                        was given. Check inputs and try again."))
                },
                "kelvin" | "k" if temperature.is_none() => {
                    temperature = Some("kelvin")
                },
                "celsius" | "c" if temperature.is_none() => {
                    temperature = Some("celsius")
                },
                "metre" if !distance => {
                    distance = true
                },
                "ev" | "joule" if energy.is_none() => {
                    energy = Some(word)
                },
                "gy" if !radiation => {
                    radiation = true
                },
                "second" | "minute" 
                | "hour" | "day" | "year" if time.is_none() => {
                    time = Some(word)
                },
                "kaannum" | "ka_annum" | "ka" | "ka_annum" if time.is_none() => {
                    time = Some("Ka")
                },
                "maannum" | "ma_annum" | "ma" | "ma_annum" if time.is_none() => {
                    time = Some("Ma")
                },
                "square" => { 
                    if power.is_none() {
                        power = Some(2)
                    } else {
                        return Err(format!("An additional power has been specified in the unit string: \n
                                    {input}. \n 
                                    Check inputs and try again."))
                    }
                },
                "cube" => { 
                    if power.is_none() {
                        power = Some(3)
                    } else {
                        return Err(format!("An additional power has been specified in the unit string: \n
                                    {input}. \n 
                                    Check inputs and try again."))
                    }
                },
                "giga" | "mega" | "kilo" | "hecto"
                | "deca" | "deci" | "centi" | "milli"
                | "micro" | "nano" | "ang" | "pico"
                | "femto" if metric.is_none() =>
                {
                    metric = Some(word);
                },
                _ => return Err(format!("An unsupported unit word: \n 
                                        {word} \n 
                                        has been entered in the unit string: \n 
                                        {input}. \n 
                                        Check inputs and try again."))             
            }
        }    
        let power= power.unwrap_or(1);
        let (metric, dimension, ) = match (time, radiation, distance, energy, temperature, power) {
            (Some(time),false, false, None, None, 1) if metric.is_none() => {
                let kind = TimeUnit::from_str(time)?;
                (
                    MetricAbbrev::STANDARD,
                    Dimension::Time {kind}
                )
            },
            (None, true, false, None, None, 1) => {
                (
                    MetricAbbrev::from_str(metric.unwrap_or(""))?,
                    Dimension::Dose
                )
            },
            (None, false, true, None, None, 1|2|3) => {
                (
                    MetricAbbrev::from_str(metric.unwrap_or(""))?,
                    match  power {
                        1 => Dimension::Length,
                        2 => Dimension::Area,
                        3 => Dimension::Volume,
                        _ => unreachable!(),
                    }
                )
            }, 
            (None, false, false, Some(energy), None, 1) => {
                (
                    MetricAbbrev::from_str(metric.unwrap_or(""))?,
                    match energy {
                        "ev" => Dimension::Energy{kind: EnergyUnit::ElectronVolt},
                        "joule" => Dimension::Energy{kind: EnergyUnit::Joule},
                        _ => unreachable!(),
                    }  
                )
            },
            (None, false, false, None, Some(temperature), 1) if metric.is_none() => {
                let kind = TemperatureUnit::from_str(temperature)?;
                (
                    MetricAbbrev::STANDARD,
                    Dimension::Temperature {kind}
                )
            },
            _ => return Err(format!("A unit has been specified that is not currently used in the program.
                                    the unit string was {input} Check inputs and try again."))              

        };
        Ok(Self {
            metric,
            dimension
        })
    }
}

impl fmt::Display for Unit {
    /// Write the unit's human-readable display name or symbol to the formatter.
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
            MetricAbbrev::ANG => {
                format!("{name}{name3}")
            },
            _ => {
                format!("{name}{name2}{name3}")
            }
        };
        
        formatter.write_str(&name)
    }
}

/// Scalar or mutable container supporting in-place affine unit conversion.
///
/// Container implementations apply the same scale and offset recursively
/// to their elements.
pub trait UnitValue {
    /// Replace each scalar value with `value * scale + offset` in place.
    ///
    /// Container implementations preserve shape and recurse through elements;
    /// `f32` results are cast back after calculation in `f64`.
    fn apply_conversion(&mut self, scale: f64, offset: f64);
}

impl UnitValue for f64 {
    /// Replace each scalar value with `value * scale + offset` in place.
    ///
    /// Container implementations preserve shape and recurse through elements;
    /// `f32` results are cast back after calculation in `f64`.
    fn apply_conversion(&mut self, scale: f64, offset: f64) {
        *self = *self * scale + offset;
    }
}
impl UnitValue for f32 {
    /// Replace each scalar value with `value * scale + offset` in place.
    ///
    /// Container implementations preserve shape and recurse through elements;
    /// `f32` results are cast back after calculation in `f64`.
    fn apply_conversion(&mut self, scale: f64, offset: f64) {
        *self = (f64::from(*self) * scale + offset) as f32;
    }
}

impl<T: UnitValue> UnitValue for Vec<T> {
    /// Replace each scalar value with `value * scale + offset` in place.
    ///
    /// Container implementations preserve shape and recurse through elements;
    /// `f32` results are cast back after calculation in `f64`.
    fn apply_conversion(&mut self, scale: f64, offset: f64) {
        for value in self {
            value.apply_conversion(scale, offset);
        }
    }
}

impl<T: UnitValue, const N: usize> UnitValue for [T; N] {
    /// Replace each scalar value with `value * scale + offset` in place.
    ///
    /// Container implementations preserve shape and recurse through elements;
    /// `f32` results are cast back after calculation in `f64`.
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
    /// Replace each scalar value with `value * scale + offset` in place.
    ///
    /// Container implementations preserve shape and recurse through elements;
    /// `f32` results are cast back after calculation in `f64`.
    fn apply_conversion(&mut self, scale: f64, offset: f64) {
        for value in self.iter_mut() {
            value.apply_conversion(scale, offset);
        }
    }
}

impl Unit {
    /// Return this dimension's stored base-unit scale or prefix power.
    ///
    /// Reciprocal dimensions use the same positive scale here; conversion logic
    /// handles whether that scale multiplies or divides the quantity.
    pub fn get_scaler(&self) -> f64 {
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
    /// Extract the time unit from a time, frequency, or dose-rate dimension.
    ///
    /// Returns an error for dimensions that contain no time unit.
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
    /// Extract the temperature unit, or return an error for a different dimension.
    pub fn get_temperature_unit(&self) -> Result<TemperatureUnit,String>{
            match self.dimension {

            Dimension::Temperature { kind }  => {
                Ok(kind)
            },
            _ => Err(format!("No TemperatureUnit in {:?}",self.dimension))

        }

    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Test to check dimensionless parsing
    fn parses_dimensionless_inputs() {
        for input in ["", " ", "\t\n", "1", " 1\t"] {
            assert_eq!(
                input.parse::<Unit>(),
                Ok(Unit {
                    metric: MetricAbbrev::STANDARD,
                    dimension: Dimension::Dimensionless,
                }),
                "input: {input:?}",
            );
        }
    }

    #[test]
    /// Test to parse simple dimensions and their aliases
    fn parses_simple_dimensions_and_aliases() {
        let cases = [
            ("metre", Dimension::Length),
            ("square metre", Dimension::Area),
            ("cube metre", Dimension::Volume),
            (
                "second",
                Dimension::Time {
                    kind: TimeUnit::Second,
                },
            ),
            (
                "minute",
                Dimension::Time {
                    kind: TimeUnit::Minute,
                },
            ),
            (
                "hour",
                Dimension::Time {
                    kind: TimeUnit::Hour,
                },
            ),
            (
                "day",
                Dimension::Time {
                    kind: TimeUnit::Day,
                },
            ),
            (
                "year",
                Dimension::Time {
                    kind: TimeUnit::Year,
                },
            ),
            (
                "ka",
                Dimension::Time {
                    kind: TimeUnit::KAnnum,
                },
            ),
            (
                "kaannum",
                Dimension::Time {
                    kind: TimeUnit::KAnnum,
                },
            ),
            (
                "ka_annum",
                Dimension::Time {
                    kind: TimeUnit::KAnnum,
                },
            ),
            (
                "ka annum",
                Dimension::Time {
                    kind: TimeUnit::KAnnum,
                },
            ),
            (
                "ma",
                Dimension::Time {
                    kind: TimeUnit::MAnnum,
                },
            ),
            (
                "maannum",
                Dimension::Time {
                    kind: TimeUnit::MAnnum,
                },
            ),
            (
                "ma_annum",
                Dimension::Time {
                    kind: TimeUnit::MAnnum,
                },
            ),
            (
                "ma annum",
                Dimension::Time {
                    kind: TimeUnit::MAnnum,
                },
            ),
            (
                "kelvin",
                Dimension::Temperature {
                    kind: TemperatureUnit::Kelvin,
                },
            ),
            (
                "k",
                Dimension::Temperature {
                    kind: TemperatureUnit::Kelvin,
                },
            ),
            (
                "celsius",
                Dimension::Temperature {
                    kind: TemperatureUnit::Celsius,
                },
            ),
            (
                "c",
                Dimension::Temperature {
                    kind: TemperatureUnit::Celsius,
                },
            ),
            (
                "joule",
                Dimension::Energy {
                    kind: EnergyUnit::Joule,
                },
            ),
            (
                "ev",
                Dimension::Energy {
                    kind: EnergyUnit::ElectronVolt,
                },
            ),
            ("gy", Dimension::Dose),
        ];

        for (input, dimension) in cases {
            // Exercise case folding and whitespace handling for every alias.
            for spelling in [
                input.to_owned(),
                format!(" \t{}\n", input.to_ascii_uppercase().replace(' ', " \t ")),
            ] {
                assert_eq!(
                    spelling.parse::<Unit>(),
                    Ok(Unit {
                        metric: MetricAbbrev::STANDARD,
                        dimension
                    }),
                    "input: {spelling:?}",
                );
            }
        }
    }

    #[test]
    /// Test to check metric prefixes and their accompanying values
    fn parses_metric_prefixes_and_length_powers() {
        let prefixes = [
            ("giga", MetricAbbrev::GIGA, 1e9_f64),
            ("mega", MetricAbbrev::MEGA, 1e6),
            ("kilo", MetricAbbrev::KILO, 1e3),
            ("hecto", MetricAbbrev::HECTO, 1e2),
            ("deca", MetricAbbrev::DECA, 1e1),
            ("deci", MetricAbbrev::DECI, 1e-1),
            ("centi", MetricAbbrev::CENTI, 1e-2),
            ("milli", MetricAbbrev::MILLI, 1e-3),
            ("micro", MetricAbbrev::MICRO, 1e-6),
            ("nano", MetricAbbrev::NANO, 1e-9),
            ("ang", MetricAbbrev::ANG, 1e-10),
            ("pico", MetricAbbrev::PICO, 1e-12),
            ("femto", MetricAbbrev::FEMTO, 1e-15),
        ];
        let powers = [
            ("", Dimension::Length, 1),
            ("square ", Dimension::Area, 2),
            ("cube ", Dimension::Volume, 3),
        ];

        for (prefix, metric, factor) in prefixes {
            for (power_word, dimension, exponent) in powers {
                let input = format!("{power_word}{prefix} metre");
                let unit = input
                    .parse::<Unit>()
                    .unwrap_or_else(|error| panic!("{input:?}: {error}"));
                assert_eq!(unit, Unit { metric, dimension }, "input: {input:?}");
                let expected = factor.powi(exponent);
                assert!(
                    (unit.get_scaler() / expected - 1.0).abs() < 1e-12,
                    "scale for {input:?}"
                );
            }
        }
    }

    #[test]
    /// Test checking reciprocal units
    fn parses_divisions_and_prefixed_energy_and_dose() {
        let cases = [
            ("1/metre", MetricAbbrev::STANDARD, Dimension::InverseLength),
            (
                "1 / square milli metre",
                MetricAbbrev::MILLI,
                Dimension::InverseArea,
            ),
            (
                "1 / cube centi metre",
                MetricAbbrev::CENTI,
                Dimension::Density,
            ),
            (
                "1/second",
                MetricAbbrev::STANDARD,
                Dimension::Frequency {
                    kind: TimeUnit::Second,
                },
            ),
            (
                "1 / ka annum",
                MetricAbbrev::STANDARD,
                Dimension::Frequency {
                    kind: TimeUnit::KAnnum,
                },
            ),
            (
                "Gy/second",
                MetricAbbrev::STANDARD,
                Dimension::DoseRate {
                    kind: TimeUnit::Second,
                },
            ),
            (
                " \tMILLI Gy\t/\nYEAR ",
                MetricAbbrev::MILLI,
                Dimension::DoseRate {
                    kind: TimeUnit::Year,
                },
            ),
            (
                "mega gy / ma annum",
                MetricAbbrev::MEGA,
                Dimension::DoseRate {
                    kind: TimeUnit::MAnnum,
                },
            ),
            ("milli gy", MetricAbbrev::MILLI, Dimension::Dose),
            (
                "kilo joule",
                MetricAbbrev::KILO,
                Dimension::Energy {
                    kind: EnergyUnit::Joule,
                },
            ),
            (
                "mega eV",
                MetricAbbrev::MEGA,
                Dimension::Energy {
                    kind: EnergyUnit::ElectronVolt,
                },
            ),
        ];

        for (input, metric, dimension) in cases {
            assert_eq!(
                input.parse::<Unit>(),
                Ok(Unit { metric, dimension }),
                "input: {input:?}"
            );
        }
    }

    #[test]
    /// Test to check bad input unit names
    fn rejects_malformed_and_unsupported_expressions() {
        let cases = [
            // Unknown, missing, or repeated components.
            "banana",
            "banana metre",
            "metre banana",
            "milli",
            "square",
            "metre metre",
            "gy gy",
            "second minute",
            "celsius kelvin",
            "joule ev",
            "milli kilo metre",
            "square square metre",
            "square cube metre",
            "1 metre",
            "1 1",
            "metre joule",
            "gy second",
            // Powers and prefixes that cannot be represented by these dimensions.
            "square second",
            "cube gy",
            "square joule",
            "cube kelvin",
            "milli second",
            "kilo year",
            "milli kelvin",
            "milli celsius",
            // Geological aliases must match the complete expression.
            "annum",
            "second annum",
            "ka annum annum",
            "annum ka",
            "ma annum second",
            // Missing operands, extra separators, and unsupported divisions.
            "/",
            "/second",
            "1/",
            "gy/ ",
            "1//second",
            "gy/second/year",
            "second/gy",
            "metre/1",
            "metre/second",
            "1/joule",
            "1/kelvin",
            "1/1",
            "gy/metre",
            "gy/gy",
            "1/milli second",
            "gy/milli second",
        ];

        for input in cases {
            assert!(
                input.parse::<Unit>().is_err(),
                "unexpectedly accepted {input:?}"
            );
        }
    }

    #[test]
    /// Test to check incorrect use values either side of slash
    fn division_errors_describe_the_problem_and_input() {
        let cases = [
            ("1//second", "Multiple '/' separators"),
            ("gy/", "Expected a unit on both sides"),
            ("/second", "Expected a unit on both sides"),
            ("second/gy", "Unsupported unit division"),
        ];

        for (input, message) in cases {
            let error = input.parse::<Unit>().unwrap_err();
            assert!(error.contains(message), "input: {input:?}, error: {error}");
            assert!(error.contains(input), "missing input context: {error}");
        }
    }

    #[test]
    /// Verify that input unit parsing distinguishes Celsius from Kelvin.
    fn temperature_quantity_units_parse_distinctly() {
        for name in ["celsius", "c", "Celsius"] {
            let unit: Unit = name.parse().unwrap();
            assert_eq!(
                unit.get_temperature_unit().unwrap(),
                TemperatureUnit::Celsius
            );
        }
        let unit: Unit = "kelvin".parse().unwrap();
        assert_eq!(
            unit.get_temperature_unit().unwrap(),
            TemperatureUnit::Kelvin
        );
    }
}
