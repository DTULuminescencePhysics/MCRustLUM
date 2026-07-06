/// This module contains constants that can be used throughout the program


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

    const GIGA_SQ: u32 = GIGA*GIGA;
    const MEGA_SQ: u32 = MEGA*MEGA;
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

    const GIGA_CUBE: u32 = GIGA_SQ*GIGA;
    const MEGA_CUBE: u32 = MEGA_SQ*MEGA;
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
    const BOLTZMANN: f32 = 1.380649e-23; /// Boltzmann constant in J/K (kg·m²/s²·K)
    const BOLTZMANN_EV: f32 = 8.617333262e-5;  /// Boltzmann constant in eV/K
    const PLANCK: f32 = 6.62607015e-34;   /// Planck constant in J·s (kg·m²/s)
    const PLANCK_BAR: f32 = 1.054571817e-34; /// Reduced Planck constant in J·s (kg·m²/s)
    const SPEED_LIGHT: i32 = 299792458;        /// Speed of light in m/s
    const ELEMENTARY_CHARGE: f32 = 1.602176634e-19;  /// Elementary charge in C
    const ELECTRON_MASS: f32 = 9.10938356e-31; /// Electron mass in kg
    const AVOGARDRO: f32 = 6.02214076e23;  /// Avogadro's number in 1/mol
    const GAS_CONSTANT: f32 = 8.314462618;  /// Gas constant in J/(mol·K)
    const VACUUM_PERM: f32 = 8.854187817e-12;  /// Vacuum permittivity in F/m 
    const EV_TO_J: f32 = 1.602176634e-19;  /// Conversion factor from eV to J
}

pub mod time{
    const SECOND: u32 = 1;
    const MINUTE: u32 = 60;
    const HOUR: u32 = 60*MINUTE;
    const DAY: u32 = 24*HOUR;
    const YEAR: u32 = DAY*365.2425;
    const K_ANNUM: u32 = metric::KILO*YEAR;
    const MA_ANNUM: u32 = metric::MEGA*YEAR;

}