//! This module contains constants that can be used throughout the program

/// This module contains conversions for SI unit prefixes
pub mod metric{
    const GIGA: u32 = 1_000_000_000;
    const MEGA: u32 = 1_000_000;
    const KILO: u32 = 1_000;
    const HECTO: u32 = 100;
    const DECA: u32 = 10;
    const DECI: f32 = 0.1;
    const CENTI: f32 = 0.01;
    const MILLI: f32 = 0.001;
    const MICRO: f32 = 0.000_001;
    const NANO: f32 = 0.000_000_001;
    const ANG: f32 = 0.0_000_000_001;
    const PICO: f32 = 0.000_000_000_001;
    const FEMTO: f32 = 0.000_000_000_000_001;

    const GIGA_SQ: u64 = GIGA as u64*GIGA as u64;
    const MEGA_SQ: u64 = MEGA as u64*MEGA as u64;
    const KILO_SQ: u32 = KILO*KILO;
    const HECTO_SQ: u32 = HECTO*HECTO;
    const DECA_SQ: u32 = DECA*DECA;
    const DECI_SQ: f32 = DECI*DECI;
    const CENTI_SQ: f32 = CENTI*CENTI;
    const MILLI_SQ: f32 = MILLI*MILLI;
    const MICRO_SQ: f32 = MICRO*MICRO;
    const NANO_SQ: f32 = NANO*NANO;
    const ANG_SQ: f32 = ANG*ANG;
    const PICO_SQ: f32 = PICO*PICO;
    const FEMTO_SQ: f32 = FEMTO*FEMTO;

    const MEGA_CUBE: u64 = MEGA_SQ*MEGA as u64;
    const KILO_CUBE: u32 = KILO_SQ*KILO;
    const HECTO_CUBE: u32 = HECTO_SQ*HECTO;
    const DECA_CUBE: u32 = DECA_SQ*DECA;
    const DECI_CUBE: f32 = DECI_SQ*DECI;
    const CENTI_CUBE: f32 = CENTI_SQ*CENTI;
    const MILLI_CUBE: f32 = MILLI_SQ*MILLI;
    const MICRO_CUBE: f32 = MICRO_SQ*MICRO;
    const NANO_CUBE: f32 = NANO_SQ*NANO;
    const ANG_CUBE: f32 = ANG_SQ*ANG;
    const PICO_CUBE: f32 = PICO_SQ*PICO;
    const FEMTO_CUBE: f32 = FEMTO_SQ*FEMTO;

}

/// This module contains useful Physical constants to be used throughout the program
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

/// Holds code for conversion between time units and ensures the correct time precision 
pub mod time{
    use std::fmt;    
    use crate::numeric::{PrecisionInput, TimeFloat, TimePrecision,};
    const SECOND: u32 = 1;
    const MINUTE: u32 = 60;
    const HOUR: u32 = 60*MINUTE;
    const DAY: u32 = 24*HOUR;
    const YEAR: u32 = 31_556_952;
    const K_ANNUM: u64 = 31_556_952_000;
    const MA_ANNUM: u64 = 31_556_952_000_000;

    #[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
    #[serde(rename_all = "snake_case")]
    pub enum TimeUnit {
        Second,
        Minute,
        Hour,
        Day,
        Year,
        KAnnum,
        MaAnnum
    }
    impl fmt::Display for TimeUnit {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            let name = match self {
                TimeUnit::Second => "second",
                TimeUnit::Minute => "minute",
                TimeUnit::Hour => "hour",
                TimeUnit::Day => "day",
                TimeUnit::Year => "year",
                TimeUnit::KAnnum => "Ka annum",
                TimeUnit::MaAnnum => "Ma annum",
            };

            formatter.write_str(name)
        }
    }
    
    #[derive(Debug, Clone, Copy)]
    pub enum TimeMultiplier {
        U32(u32),
        U64(u64),
    }

    impl TimeMultiplier {
        pub fn get_float_precision(self) -> TimeFloat {
            match self {
                TimeMultiplier::U32(value) => value as TimeFloat,
                TimeMultiplier::U64(value) => value as TimeFloat,
            }
        }
    }

    pub fn unit_multiplier(unit: TimeUnit) -> Option<TimeMultiplier> {
        match unit {
            TimeUnit::Second  => Some(TimeMultiplier::U32(SECOND)),
            TimeUnit::Minute  => Some(TimeMultiplier::U32(MINUTE)),
            TimeUnit::Hour  => Some(TimeMultiplier::U32(HOUR)),
            TimeUnit::Day  => Some(TimeMultiplier::U32(DAY)),
            TimeUnit::Year  => Some(TimeMultiplier::U32(YEAR)),
            TimeUnit::KAnnum  => Some(TimeMultiplier::U64(K_ANNUM)),
            TimeUnit::MaAnnum  => Some(TimeMultiplier::U64(MA_ANNUM)),
        }
       
    }
    pub fn is_ka_or_ma(unit: TimeUnit) -> bool {
        match unit {
            TimeUnit::KAnnum  => true,
            TimeUnit::MaAnnum  => true,
            _ => false,
        }
    }
    pub fn convert_to_seconds<T>(unit: TimeUnit, value: T) -> Option<T::Output>
    where
        T: PrecisionInput<TimePrecision>,
    {
        Some(value.multiply_to_precision(
       unit_multiplier(unit)?
                  .get_float_precision()
            )
        )
    }
}
/// Contains Temperature conversion behaviour 
pub mod temperature {
    use crate::numeric::{ElementWise, Float, PrecisionInput, ProgramPrecision,};
    #[derive(Debug, Clone, Copy, PartialEq, serde::Deserialize)]
    #[serde(rename_all = "snake_case")]
    pub enum TemperatureUnit {
        Celsius,
        Kelvin,
    }

    /// Convert a scalar or container of temperatures to kelvin at program precision.
    pub fn convert_to_kelvin<T>(unit: TemperatureUnit, value: T) -> Option<T::Output>
    where
        T: PrecisionInput<ProgramPrecision>,
        T::Output: ElementWise<Float, Output = T::Output>,
    {
        const CELSIUS_TO_KELVIN: Float = 273.15;

        let value = value.map_to_precision();

        match unit {
            TemperatureUnit::Celsius => value.element_add(&CELSIUS_TO_KELVIN),
            TemperatureUnit::Kelvin => Some(value),
        }
    }

        /// Convert a scalar or container of temperatures to kelvin at program precision.
    pub fn convert_to_celsius<T>(unit: TemperatureUnit, value: T) -> Option<T::Output>
    where
        T: PrecisionInput<ProgramPrecision>,
        T::Output: ElementWise<Float, Output = T::Output>,
    {
        const CELSIUS_TO_KELVIN: Float = 273.15;

        let value = value.map_to_precision();

        match unit {
            TemperatureUnit::Kelvin  => value.element_sub(&CELSIUS_TO_KELVIN),
            TemperatureUnit::Celsius => Some(value),
        }
    }

}

#[cfg(test)]
mod tests {
    use super::temperature::{TemperatureUnit, convert_to_kelvin};
    use crate::numeric::Float;

    #[test]
    fn converts_celsius_to_kelvin_at_program_precision() {
        let temperature = convert_to_kelvin(TemperatureUnit::Celsius, 20.0_f32).unwrap();

        assert_eq!(temperature, 293.15 as Float);
    }

    #[test]
    fn converts_temperature_vectors_to_kelvin() {
        let temperatures =
            convert_to_kelvin(TemperatureUnit::Celsius, vec![0_i32, 100_i32]).unwrap();

        assert_eq!(temperatures, vec![273.15 as Float, 373.15 as Float]);
    }

    #[test]
    fn preserves_kelvin_values_while_converting_precision() {
        let temperature = convert_to_kelvin(TemperatureUnit::Kelvin, 300.0_f32).unwrap();

        assert_eq!(temperature, 300.0 as Float);
    }
}


