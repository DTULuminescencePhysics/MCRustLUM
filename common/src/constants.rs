// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! SI conversions and physical constants used by the transport equations.
/// Internal scalar, square, and cubic SI prefix factors.
///
/// These factors are used by unit parsing and quantity conversion.
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

    /// Cubed giga factor for volume conversions.
    const GIGA_CUBE: f64 = GIGA_SQ*GIGA;
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

    /// Metric prefix used to scale physical quantities, including the ångström length scale.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum MetricAbbrev {
        /// Giga prefix (`10^9`).
        GIGA,
        /// Mega prefix (`10^6`).
        MEGA,
        /// Kilo prefix (`10^3`).
        KILO,
        /// Hecto prefix (`10^2`).
        HECTO,
        /// Deca prefix (`10^1`).
        DECA,
        /// Deci prefix (`10^-1`).
        DECI,
        /// Centi prefix (`10^-2`).
        CENTI,
        /// Base unit without a prefix.
        STANDARD,
        /// Milli prefix (`10^-3`).
        MILLI,
        /// Micro prefix (`10^-6`).
        MICRO,
        /// Nano prefix (`10^-9`).
        NANO,
        /// Ångström length scale (`10^-10` metres).
        ANG,
        /// Pico prefix (`10^-12`).
        PICO,
        /// Femto prefix (`10^-15`).
        FEMTO,
    }
     impl FromStr for MetricAbbrev {
        type Err = String;

        /// Parse a case-insensitive prefix name; unrecognized names select `STANDARD`.
        ///
        /// This parser recognizes the names explicitly listed in its match table;
        /// it does not reject unknown text.
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s.to_lowercase().as_str() {
                "giga"   => Ok(Self::GIGA),
                "mega"   => Ok(Self::MEGA),
                "kilo"   => Ok(Self::KILO),
                "hecto"  => Ok(Self::HECTO),
                "deca"  => Ok(Self::DECA),
                "deci" => Ok(Self::DECI),
                "centi" => Ok(Self::CENTI),
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


    impl MetricAbbrev {

        /// Return the display symbol for this prefix, or an empty string for the base unit.
        pub fn get_name(&self)-> &'static str {
            match self {
                MetricAbbrev::STANDARD => "",
                MetricAbbrev::GIGA => "G",
                MetricAbbrev::MEGA => "M",
                MetricAbbrev::KILO => "k",
                MetricAbbrev::HECTO => "h",
                MetricAbbrev::DECA => "da",
                MetricAbbrev::DECI => "d",
                MetricAbbrev::CENTI => "c",
                MetricAbbrev::MILLI => "m",
                MetricAbbrev::MICRO => "μ",
                MetricAbbrev::NANO => "n",
                MetricAbbrev::ANG => "Å",
                MetricAbbrev::PICO => "p",
                MetricAbbrev::FEMTO => "f",
            }
        }
        /// Return the configured linear prefix multiplier relative to the base unit.
        pub fn get_metric(&self)-> f64 {
            match self {
                MetricAbbrev::STANDARD => 1.0,
                MetricAbbrev::GIGA => GIGA,
                MetricAbbrev::MEGA => MEGA,
                MetricAbbrev::KILO => KILO,
                MetricAbbrev::HECTO => HECTO,
                MetricAbbrev::DECA => DECA,
                MetricAbbrev::DECI => DECI,
                MetricAbbrev::CENTI => CENTI,
                MetricAbbrev::MILLI => MILLI,
                MetricAbbrev::MICRO => MICRO,
                MetricAbbrev::NANO => NANO,
                MetricAbbrev::ANG => ANG,
                MetricAbbrev::PICO => PICO,
                MetricAbbrev::FEMTO => FEMTO,
            }
        }
        /// Return the configured squared-prefix multiplier used by area conversions.
        pub fn get_square_metric(&self)-> f64 {
            match self {
                MetricAbbrev::STANDARD => 1.0,
                MetricAbbrev::GIGA => GIGA_SQ,
                MetricAbbrev::MEGA => MEGA_SQ,
                MetricAbbrev::KILO => KILO_SQ,
                MetricAbbrev::HECTO => HECTO_SQ,
                MetricAbbrev::DECA => DECA_SQ,
                MetricAbbrev::DECI => DECI_SQ,
                MetricAbbrev::CENTI => CENTI_SQ,
                MetricAbbrev::MILLI => MILLI_SQ,
                MetricAbbrev::MICRO => MICRO_SQ,
                MetricAbbrev::NANO => NANO_SQ,
                MetricAbbrev::ANG => ANG_SQ,
                MetricAbbrev::PICO => PICO_SQ,
                MetricAbbrev::FEMTO => FEMTO_SQ,
            }
        }
        /// Return the cubed prefix multiplier for volume conversions.
        pub fn get_cube_metric(&self)-> f64 {
            match self {
                MetricAbbrev::STANDARD => 1.0,
                MetricAbbrev::GIGA => GIGA_CUBE,
                MetricAbbrev::MEGA => MEGA_CUBE,
                MetricAbbrev::KILO => KILO_CUBE,
                MetricAbbrev::HECTO => HECTO_CUBE,
                MetricAbbrev::DECA => DECA_CUBE,
                MetricAbbrev::DECI => DECI_CUBE,
                MetricAbbrev::CENTI => CENTI_CUBE,
                MetricAbbrev::MILLI => MILLI_CUBE,
                MetricAbbrev::MICRO => MICRO_CUBE,
                MetricAbbrev::NANO => NANO_CUBE,
                MetricAbbrev::ANG => ANG_CUBE,
                MetricAbbrev::PICO => PICO_CUBE,
                MetricAbbrev::FEMTO => FEMTO_CUBE,
            }
        }

    }

}

/// Fundamental constants expressed in SI units, plus Boltzmann's constant in eV/K.
pub mod physical_constants{
    use crate::numeric::{ElementWise, Float, PrecisionInput, ProgramPrecision,};

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

    /// Energy unit attached to an energy dimension.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum EnergyUnit {
        /// Energy in joules.
        Joule,
        /// Energy in electron volts.
        ElectronVolt,
    }
    impl EnergyUnit {
        /// Return `J` or `eV` for energy labels.
        pub fn get_name(&self,) -> &'static str {
            match self {
                EnergyUnit::Joule => "J",
                EnergyUnit::ElectronVolt => "eV",
               
            }
        }
        /// Convert values from energy unit `self` to `unit`, preserving container shape.
        ///
        /// Electron volts to joules multiply by `EV_TO_J`; the reverse divides by it.
        /// Equal units retain values. Returns the optional element-wise result.
        pub fn convert_unit<T>(&self, unit: EnergyUnit, value:T)-> Option<T::Output>
        where
            T: PrecisionInput<ProgramPrecision>,
            T::Output: ElementWise<Float, Output = T::Output>,
        {
            let value = value.map_to_precision();
            match (self, unit) {     
                (EnergyUnit::ElectronVolt, EnergyUnit::Joule) => {
                    value.element_mul(&EV_TO_J)
                },
                (EnergyUnit::Joule, EnergyUnit::ElectronVolt) => {
                    value.element_div(&EV_TO_J)
                },
                (EnergyUnit::Joule, EnergyUnit::Joule) => Some(value),
                (EnergyUnit::ElectronVolt, EnergyUnit::ElectronVolt) =>  Some(value),
            } 
        }
    }

}

/// Time-unit conversion at the precision used by the simulation clock.
pub mod time{
    use std::{fmt::{self, Debug}, str::FromStr};    
    use crate::numeric::{ElementWise, PrecisionInput, TimeFloat, TimePrecision,};
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
    const M_ANNUM: f64 = 31_556_952_000_000.0;

    /// Unit used for time values supplied to the simulation.
    ///
    /// Geological ages use annum-based names: `KAnnum` is one thousand years
    /// and `MAnnum` is one million years.
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
        MAnnum
    }
    impl FromStr for TimeUnit {
        type Err = String;

        /// Parse a time-unit name or supported geological-age alias; reject unknown names.
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s.to_lowercase().as_str() {
                "second"   => Ok(Self::Second),
                "minute"   => Ok(Self::Minute),
                "hour"     => Ok(Self::Hour),
                "day"      => Ok(Self::Day),
                "year"     => Ok(Self::Year),
                "ka annum" | "kaannum" | "ka_annum" | "ka" | "kannum" => Ok(Self::KAnnum),
                "ma annum" | "maannum" | "ma_annum" | "ma" | "mannum" => Ok(Self::MAnnum),
                _ => Err(format!("Invalid TimeUnit: {s}")),
            }
        }
    }
    impl fmt::Display for TimeUnit {
        /// Write the unit's human-readable display name or symbol to the formatter.
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            let name = match self {
                TimeUnit::Second => "second",
                TimeUnit::Minute => "minute",
                TimeUnit::Hour => "hour",
                TimeUnit::Day => "day",
                TimeUnit::Year => "year",
                TimeUnit::KAnnum => "Ka",
                TimeUnit::MAnnum => "Ma",
            };

            formatter.write_str(name)
        }
    }
    impl TimeUnit{
        /// Checks if the time unit is Ka or Ma and so ordering is reversed
        pub fn is_ka_or_ma(&self) -> bool {
            match self {
                TimeUnit::KAnnum  => true,
                TimeUnit::MAnnum  => true,
                _ => false,
            }
        }
        /// Return the number of seconds in one `unit`.
        pub fn get_time_unit_metric(&self,) -> TimeFloat {
            match self {
                TimeUnit::Second => SECOND as TimeFloat,
                TimeUnit::Minute => MINUTE as TimeFloat,
                TimeUnit::Hour => HOUR as TimeFloat,
                TimeUnit::Day => DAY as TimeFloat,
                TimeUnit::Year => YEAR as TimeFloat,
                TimeUnit::KAnnum => K_ANNUM as TimeFloat,
                TimeUnit::MAnnum => M_ANNUM as TimeFloat,
            }
        }
        /// Return the human-readable time-unit name used in labels.
        pub fn get_name(&self) -> &'static str {
            match self {
                TimeUnit::Second => "second",
                TimeUnit::Minute => "minute",
                TimeUnit::Hour => "hour",
                TimeUnit::Day => "day",
                TimeUnit::Year => "year",
                TimeUnit::KAnnum => "Ka",
                TimeUnit::MAnnum => "Ma",
            }
        }
        /// Convert a scalar or container from `unit` into the target time unit `self`.
        ///
        /// Scales by source-seconds divided by target-seconds at time precision,
        /// preserving container shape. Returns the element-wise operation's optional result.
        pub fn convert_unit<T>(&self, unit: &TimeUnit, value:T)-> Option<T::Output>
        where
            T: PrecisionInput<TimePrecision>,
            T::Output: ElementWise<TimeFloat, Output = T::Output>,
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
        /// Write the unit's human-readable display name or symbol to the formatter.
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

        /// Parse Kelvin/Celsius names or symbols without case sensitivity; reject unknown units.
        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s.to_lowercase().as_str() {
                "kelvin" | "k" => Ok(Self::Kelvin),
                "celsius" | "°c" => Ok(Self::Celsius),
                _ => Err(format!("Invalid TemperatureUnit: {s}")),
            }
        }
    }
    impl TemperatureUnit {
        /// Convert values using the source/target pair `(self, unit)` at program precision.
        ///
        /// The Celsius/Kelvin pair adds 273.15 and the reverse pair subtracts it;
        /// equal units retain values. Container shape is preserved and conversion
        /// returns the optional element-wise result.
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
        /// Return `°C` or `K` for temperature labels.
        pub fn get_name(&self) -> &'static str {
            match self {
                TemperatureUnit::Celsius => "°C",
                TemperatureUnit::Kelvin => "K",
            }
        }
    }

}




#[cfg(test)]
mod tests {
    use super::temperature::TemperatureUnit;
    use super::time::TimeUnit;
    use crate::numeric::Float;

   

    #[test]
    /// Verify that converts celsius to kelvin at program precision.
    fn converts_celsius_to_kelvin_at_program_precision() {
        let temperature = TemperatureUnit::Celsius.convert_unit(&TemperatureUnit::Kelvin, 20.0_f32).unwrap();

        assert_eq!(temperature, 293.15 as Float);
    }

    #[test]
    /// Verify that converts temperature vectors to kelvin.
    fn converts_temperature_vectors_to_kelvin() {
        let temperatures = TemperatureUnit::Celsius.convert_unit(&TemperatureUnit::Kelvin, vec![0_i32, 100_i32]).unwrap();


        assert_eq!(temperatures, vec![273.15 as Float, 373.15 as Float]);
    }

    #[test]
    /// Verify that preserves kelvin values while converting precision.
    fn preserves_kelvin_values_while_converting_precision() {
        let temperature = TemperatureUnit::Kelvin.convert_unit(&TemperatureUnit::Kelvin, 300.0_f32).unwrap();


        assert_eq!(temperature, 300.0 as Float);
    }

    #[test]
    /// Verify that converts kelvin to celsius for scalars and vectors.
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
    /// Verify that temperature units parse and display consistently.
    fn temperature_units_parse_and_display_consistently() {
        assert_eq!("kelvin".parse(), Ok(TemperatureUnit::Kelvin));
        assert_eq!("celsius".parse(), Ok(TemperatureUnit::Celsius));
        assert_eq!(TemperatureUnit::Kelvin.to_string(), "K");
        assert_eq!(TemperatureUnit::Celsius.to_string(), "°C");
    }

    #[test]
    /// Verify that time units display with human readable names.
    fn time_units_display_with_human_readable_names() {
        let names = [
            (TimeUnit::Second, "second"),
            (TimeUnit::Minute, "minute"),
            (TimeUnit::Hour, "hour"),
            (TimeUnit::Day, "day"),
            (TimeUnit::Year, "year"),
            (TimeUnit::KAnnum, "Ka"),
            (TimeUnit::MAnnum, "Ma"),
        ];

        for (unit, expected) in names {
            assert_eq!(unit.to_string(), expected);
        }
    }

    #[test]
    /// Verify that time multipliers and seconds cover short and geological units.
    fn time_multipliers_and_seconds_cover_short_and_geological_units() {
        assert!(matches!(
            TimeUnit::Hour.get_time_unit_metric(),
            3_600.0,
        ));
        assert!(matches!(
            TimeUnit::KAnnum.get_time_unit_metric(),
            31_556_952_000.0,
        ));
        assert!(!(TimeUnit::Year.is_ka_or_ma()));
        assert!((TimeUnit::KAnnum.is_ka_or_ma()));
        assert!((TimeUnit::MAnnum.is_ka_or_ma()));

        assert_eq!(
            TimeUnit::Second.convert_unit(&TimeUnit::Day,  vec![0_i32, 1, 2]),
            Some(vec![0.0, 86_400.0, 172_800.0]),
        );
        assert_eq!(
            TimeUnit::Second.convert_unit(&TimeUnit::MAnnum,  1.0_f64),
            Some(31_556_952_000_000.0),
        );
    }
}
