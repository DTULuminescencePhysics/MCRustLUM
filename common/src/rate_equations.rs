// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only 

//! Generic rate equations used by the Monte Carlo luminescence model.
//!
//! All functions operate through the element-wise traits in `numeric.rs`,
//! allowing the same equation to accept scalars, vectors, or ndarrays.
//! Mixed input shapes are supported when a corresponding broadcasting
//! implementation exists, such as vector energies with scalar temperature.
//!
//! Most public functions return `Option<T>`. `None` indicates that an
//! element-wise operation could not be completed, usually because input
//! container shapes were incompatible. Successful rate values are converted
//! to `TimePrecision` before being returned.

use crate::numeric::{Float, ElementWise, ElementWiseUnary, PrecisionInput, TimePrecision};
use crate::constants::physical_constants::{BOLTZMANN_EV};


/// Function to calculate the exponential of the energy over KbT.
///
/// Energy and temperature may use different container types. For example,
/// `e` may be a `Vec<Float>` while `temp` is a single `Float`, provided the
/// required `ElementWise` implementations exist.
pub fn exponential_energy_over_kb_t<E, T, ENeg, KT, Ratio, Exp>(
    e: &E,
    temp: &T,
) -> Option<Exp>
where
    E: ElementWise<Float, Output = ENeg>,
    T: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
{
    let negative_energy = e.element_mul(&-1.0)?;
    let kb_t = temp.element_mul(&BOLTZMANN_EV)?;
    Some(negative_energy.element_div(&kb_t)?.element_exp())
}

/// Calculate the first-order delocalised release rate.
/// The model evaluates the per-particle rate `s * exp(-E / (k_B T))`, where
/// `s` is the attempt frequency. Population scaling is applied by callers that
/// calculate population changes.
pub fn first_order_delocalised_rate_equation<E, S, T, ENeg, KT, Ratio, Exp, V>(
    e_cb: &E,
    s_frequency: &S,
    temp: &T,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    E: ElementWise<Float, Output = ENeg>,
    T: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<S, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    let exponent = exponential_energy_over_kb_t(e_cb, temp)?;
    Some(exponent.element_mul(s_frequency)?.map_to_precision())
}

/// Calculate the second-order delocalised release rate.
///
/// The model evaluates the per-particle rate `s * exp(-E / (k_B T))`.
pub fn second_order_delocalised_rate_equation<E, S, T, ENeg, KT, Ratio, Exp, V>(
    e_cb: &E,
    s_frequency: &S,
    temp: &T,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    E: ElementWise<Float, Output = ENeg>,
    T: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<S, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    let exponent = exponential_energy_over_kb_t(e_cb, temp)?;
    Some(exponent.element_mul(s_frequency)?.map_to_precision())
}

/// Calculate a general-order delocalised release rate.
///
/// The model evaluates the per-particle rate `s * exp(-E / (k_B T))`.
/// `order` remains part of the configured model identity but population-order
/// scaling is no longer performed by this function.
pub fn general_order_delocalised_rate_equation<E, S, Temp, Order, ENeg, KT, Ratio, Exp, V>(
    e_cb: &E,
    s_frequency: &S,
    temp: &Temp,
    _order: &Order,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<S, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    let exponent = exponential_energy_over_kb_t(e_cb, temp)?;
    Some(exponent.element_mul(s_frequency)?.map_to_precision())
}

/// Calculate recombination loss from the hole population.
///
/// The returned magnitude is `nc * m * recomb`, where `nc` is the
/// conduction-band population, `m` is the occupied-hole population, and
/// `recomb` is the recombination coefficient.
pub fn hole_change_delocalised_rate_equation<M, Nc, Recomb, MNc, V>(
    m: &M,
    nc: &Nc,
    recomb: &Recomb,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    Nc: ElementWise<M, Output = MNc>,
    MNc: ElementWise<Recomb, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    Some(
        nc.element_mul(m)?
            .element_mul(recomb)?
            .map_to_precision(),
    )
}

/// Calculate conduction-band retrapping into available traps.
///
/// Available traps are `n_tot - n`, giving the rate
/// `nc * (n_tot - n) * retrap`.
pub fn retrapping_change_delocalised_rate_equation<NTot, N, Nc, Retrap, Available, NcAvailable, V>(
    n_tot: &NTot,
    n: &N,
    nc: &Nc,
    retrap: &Retrap,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    NTot: ElementWise<N, Output = Available>,
    Nc: ElementWise<Available, Output = NcAvailable>,
    NcAvailable: ElementWise<Retrap, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    let available = n_tot.element_sub(n)?;
    Some(
        nc.element_mul(&available)?
            .element_mul(retrap)?
            .map_to_precision(),
    )
}

/// Calculate the net change in occupied trap concentration.
///
/// Thermal release decreases the trap population, while conduction-band
/// retrapping increases it. The returned value is
/// `retrapping_rate - thermal_release_rate`.
pub fn trap_change_delocalised_rate_equation
< E, S, N, Temp, NTot, Nc, Retrap, ENeg, KT, Ratio, Exp, ThermalRaw,
  ThermalBase, PopulationRaw, PopulationOut, ThermalOut, Available, NcAvailable,
  RetrapRaw, RetrapOut, V, >(
    e_cb: &E,
    s_frequency: &S,
    n: &N,
    temp: &Temp,
    n_tot: &NTot,
    nc: &Nc,
    retrap: &Retrap,
) -> Option<V >
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<S, Output = ThermalRaw>,
    ThermalRaw: PrecisionInput<TimePrecision, Output = ThermalBase>,
    N: ElementWise<Float, Output = PopulationRaw>,
    PopulationRaw: PrecisionInput<TimePrecision, Output = PopulationOut>,
    ThermalBase: ElementWise<PopulationOut, Output = ThermalOut>,
    NTot: ElementWise<N, Output = Available>,
    Nc: ElementWise<Available, Output = NcAvailable>,
    NcAvailable: ElementWise<Retrap, Output = RetrapRaw>,
    RetrapRaw: PrecisionInput<TimePrecision, Output = RetrapOut>,
    RetrapOut: ElementWise<ThermalOut, Output = V>,

{
    let per_particle_to_cb = first_order_delocalised_rate_equation(
        e_cb,
        s_frequency,
        temp,
    )?;
    let population = n.element_mul(&1.0)?.map_to_precision();
    let to_cb = per_particle_to_cb.element_mul(&population)?;
    let back_cb = retrapping_change_delocalised_rate_equation(
        n_tot,
        n,
        nc,
        retrap,
    )?;

    back_cb.element_sub(&to_cb)
}

// /// Calculate the net change in conduction-band carrier concentration.
///
/// Thermal release supplies carriers to the conduction band, while
/// retrapping and recombination remove them. This implementation follows
/// the existing simulation sign convention and returns
/// `retrapping - recombination - thermal_release`. 
pub fn cb_band_change_delocalised_rate_equation
< E, S, N, Temp, NTot, Nc, M, Retrap, Recomb, ENeg, KT, Ratio, Exp, 
  ThermalRaw, ThermalBase, PopulationRaw, PopulationOut, ThermalOut, Available,
  NcAvailable, RetrapRaw, RetrapOut, MNc, HoleRaw, HoleOut, LossesRaw, Losses, V >(
    e_cb: &E,
    s_frequency: &S,
    n: &N,
    temp: &Temp,
    n_tot: &NTot,
    nc: &Nc,
    m: &M,
    retrap: &Retrap,
    recomb: &Recomb,
) -> Option< V >
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<S, Output = ThermalRaw>,
    ThermalRaw: PrecisionInput<TimePrecision, Output = ThermalBase>,
    N: ElementWise<Float, Output = PopulationRaw>,
    PopulationRaw: PrecisionInput<TimePrecision, Output = PopulationOut>,
    ThermalBase: ElementWise<PopulationOut, Output = ThermalOut>,
    NTot: ElementWise<N, Output = Available>,
    Nc: ElementWise<Available, Output = NcAvailable>
        + ElementWise<M, Output = MNc>,
    NcAvailable: ElementWise<Retrap, Output = RetrapRaw>,
    RetrapRaw: PrecisionInput<TimePrecision, Output = RetrapOut>,
    MNc: ElementWise<Recomb, Output = HoleRaw>,
    HoleRaw: PrecisionInput<TimePrecision, Output = HoleOut>,
    HoleOut: ElementWise<ThermalOut, Output = LossesRaw>,
    LossesRaw: PrecisionInput<TimePrecision, Output = Losses>,
    RetrapOut: ElementWise<Losses, Output = V>,
    

{
    let per_particle_from_trap = first_order_delocalised_rate_equation(
        e_cb,
        s_frequency,
        temp,
    )?;
    let population = n.element_mul(&1.0)?.map_to_precision();
    let from_trap = per_particle_from_trap.element_mul(&population)?;
    let to_trap = retrapping_change_delocalised_rate_equation(
        n_tot,
        n,
        nc,
        retrap,
    )?;
    let to_hole = hole_change_delocalised_rate_equation(
        m,
        nc,
        recomb,
    )?;

    let losses = to_hole.element_add(&from_trap)?.map_to_precision();
    to_trap.element_sub(&losses)
}

/// Calculate the quasi-equilibrium delocalised recombination rate.
///
/// The thermal release rate is multiplied by the fraction of released
/// carriers expected to recombine rather than retrap.
pub fn quasi_equ_delocalised_rate_equation
< E, S, N, Temp, NTot, M, Retrap, Recomb, ENeg, KT, Ratio, Exp, ThermalRaw,
  ThermalBase, PopulationRaw, PopulationOut, ThermalOut, Available, ToTrap,
  ToHole, Total, FractionRaw, FractionOut, V >(
    e_cb: &E,
    s_frequency: &S,
    n: &N,
    temp: &Temp,
    n_tot: &NTot,
    m: &M,
    retrap: &Retrap,
    recomb: &Recomb,
) -> Option< V >
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<S, Output = ThermalRaw>,
    ThermalRaw: PrecisionInput<TimePrecision, Output = ThermalBase>,
    N: ElementWise<Float, Output = PopulationRaw>,
    PopulationRaw: PrecisionInput<TimePrecision, Output = PopulationOut>,
    ThermalBase: ElementWise<PopulationOut, Output = ThermalOut>,
    NTot: ElementWise<N, Output = Available>,
    Retrap: ElementWise<Available, Output = ToTrap>,
    M: ElementWise<Recomb, Output = ToHole>,
    ToTrap: ElementWise<ToHole, Output = Total>,
    ToHole: ElementWise<Total, Output = FractionRaw>,
    FractionRaw: PrecisionInput<TimePrecision, Output = FractionOut>,
    ThermalOut: ElementWise<FractionOut, Output = V>,
{
    let per_particle_to_cb = first_order_delocalised_rate_equation(
        e_cb,
        s_frequency,
        temp,
    )?;
    let population = n.element_mul(&1.0)?.map_to_precision();
    let to_cb = per_particle_to_cb.element_mul(&population)?;
    let available = n_tot.element_sub(n)?;
    let to_trap = retrap.element_mul(&available)?;
    let to_hole = m.element_mul(recomb)?;
    let total = to_trap.element_add(&to_hole)?;
    let fraction = to_hole.element_div(&total)?.map_to_precision();

    to_cb.element_mul(&fraction)
}


/// Localised Transitions 
/// Calculate a localised tunnelling rate with exponential distance decay.
///
/// The model evaluates `b * exp(-alpha * r)`, where `alpha` is the decay
/// constant, `b` is the attempt frequency, and `r` is the tunnelling distance.
pub fn tunnelling_rate<Alpha, B, R, NegativeAlpha, Product, Exp, V>(
    alpha: &Alpha,
    b: &B,
    r: &R,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    Alpha: ElementWise<Float, Output = NegativeAlpha>,
    NegativeAlpha: ElementWise<R, Output = Product>,
    Product: ElementWiseUnary<Output = Exp>,
    B: ElementWise<Exp, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    let negative_alpha = alpha.element_mul(&-1.0)?;
    let exponent = negative_alpha.element_mul(r)?.element_exp();
    Some(b.element_mul(&exponent)?.map_to_precision())
}

/// Filling 
/// Calculate filling of currently available traps under an external dose.
///
/// The model evaluates `(d_dot / d0) * (n_tot - n)`, where `d0` is the
/// characteristic dose and `d_dot` is the applied dose rate.
pub fn filling_rate<D0, DDot, N, NTot, DoseRatio, Available, V>(
    d0: &D0,
    d_dot: &DDot,
    n: &N,
    n_tot: &NTot,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    DDot: ElementWise<D0, Output = DoseRatio>,
    NTot: ElementWise<N, Output = Available>,
    DoseRatio: ElementWise<Available, Output = V>,
    V: PrecisionInput<TimePrecision>,
{
    let dose_ratio = d_dot.element_div(d0)?;
    let available = n_tot.element_sub(n)?;
    Some(dose_ratio.element_mul(&available)?.map_to_precision())
}

/// Calculate trap filling when newly released carriers may recombine directly.
///
/// The basic filling rate is multiplied by the affinity-derived probability
/// that the carrier is retrapped rather than recombined.   
pub fn filling_with_recombination_rate
< D0, DDot, N, NTot, M, Retrap, Recomb, DoseRatio, Available, FillRaw, 
  FillOut, ToTrap, ToHole, Total,FractionRaw, Fraction, V >(
    d0: &D0,
    d_dot: &DDot,
    n: &N,
    n_tot: &NTot,
    m: &M,
    retrap: &Retrap,
    recomb: &Recomb,
) -> Option<V >
where
    DDot: ElementWise<D0, Output = DoseRatio>,
    NTot: ElementWise<N, Output = Available>,
    DoseRatio: ElementWise<Available, Output = FillRaw>,
    FillRaw: PrecisionInput<TimePrecision, Output = FillOut>,
    Retrap: ElementWise<Available, Output = ToTrap>,
    M: ElementWise<Recomb, Output = ToHole>,
    ToTrap: ElementWise<ToHole, Output = Total>
        + ElementWise<Total, Output = FractionRaw>,
    FractionRaw: PrecisionInput<TimePrecision, Output = Fraction>,
    FillOut: ElementWise<Fraction, Output = V>,
   
{
    let fill = filling_rate(d0, d_dot, n, n_tot)?;
    let fraction = retrapping_by_affinity(n, n_tot, m, retrap, recomb)?;
    fill.element_mul(&fraction)
}

/// Return the affinity-derived probability of retrapping.
///
/// The result is `to_trap / (to_trap + to_hole)`.
pub fn retrapping_by_affinity 
< N, NTot, M, Retrap, Recomb, Available, 
  ToTrap, ToHole, Total,FractionRaw, W > (
    n: &N,
    n_tot: &NTot,
    m: &M,
    retrap: &Retrap,
    recomb: &Recomb,
) -> Option < W >
where
    NTot: ElementWise<N, Output = Available>,
    Retrap: ElementWise<Available, Output = ToTrap>,
    M: ElementWise<Recomb, Output = ToHole>,
    ToTrap: ElementWise<ToHole, Output = Total>
        + ElementWise<Total, Output = FractionRaw>,
    FractionRaw: PrecisionInput<TimePrecision, Output = W>,
{   
    let available = n_tot.element_sub(n)?;
    let to_trap = retrap.element_mul(&available)?;
    let to_hole = m.element_mul(recomb)?;
    let total = to_trap.element_add(&to_hole)?;
    Some(to_trap.element_div(&total)?.map_to_precision())

}
/// Return the affinity-derived probability of recombination.
///
/// The result is `to_hole / (to_trap + to_hole)`. It complements
/// `retrapping_by_affinity` when the denominator is nonzero.
pub fn recombination_by_affinity 
< N, NTot, M, Retrap, Recomb, Available, 
  ToTrap, ToHole, Total,FractionRaw, W > (
    n: &N,
    n_tot: &NTot,
    m: &M,
    retrap: &Retrap,
    recomb: &Recomb,
) -> Option < W >
where
    NTot: ElementWise<N, Output = Available>,
    Retrap: ElementWise<Available, Output = ToTrap>,
    M: ElementWise<Recomb, Output = ToHole>,
    ToHole: ElementWise<ToTrap, Output = Total>
        + ElementWise<Total, Output = FractionRaw>,
    FractionRaw: PrecisionInput<TimePrecision, Output = W>,
{   
    let available = n_tot.element_sub(n)?;
    let to_trap = retrap.element_mul(&available)?;
    let to_hole = m.element_mul(recomb)?;
    let total =to_hole.element_add(&to_trap)?;
    Some(to_hole.element_div(&total)?.map_to_precision())

}
/// Calculate the current distance-dependent retrapping factor.
///
/// The implementation evaluates `out = prefactor * exp((r / mu)^2)`. The name is
/// retained for compatibility, although the function returns the model
/// reciprocal factor to be used to select a transition tau = -ln(u) * out.
pub fn retrapping_probability_by_r <Pf, Mu, R, Frac, Sq, Exp, V> (
    prefactor: &Pf,
    mu: &Mu,
    r: &R 
) ->  Option<<V as PrecisionInput<TimePrecision>>::Output>
where 
    R: ElementWise<Mu, Output = Frac>,
    Frac: ElementWiseUnary<Output = Sq >,
    Sq: ElementWiseUnary< Output = Exp >,
    Pf: ElementWise<Exp, Output = V >,
   
    V: PrecisionInput<TimePrecision>,

{
    let frac = r.element_div(mu)?; 
    let sq = frac.element_powf(2);
    let e = sq.element_exp();
    
    Some(prefactor.element_mul(&e)?.map_to_precision())
}

/// Internal transitions
/// Calculate the thermal occupation weights of the ground and excited states.
///
/// For `x = s_e*exp(-E / (k_B T))`, the returned pair is
/// `(ground_weight, excited_weight) = (s_g / (s_g + x), x / (s_g + x))`.
/// Consequently, the two weights sum to one for every element. Energy and
/// temperature accept the same scalar, vector, and array combinations as the
/// thermal state-rate functions below.
pub fn ground_excited_state_weights
<E, SE, SG, Temp, ENeg, KT, Ratio, Exp, ExcitationRaw,
  ExcitationBase, Denominator, SGM, SGOut,      

GroundWeight, ExcitedWeight>(
    e: &E,
    s_frequency_e: &SE,
    s_frequency_g: &SG,
    temp: &Temp,
) -> Option<(GroundWeight, ExcitedWeight)>
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<SE, Output = ExcitationRaw>,
    ExcitationRaw: PrecisionInput<TimePrecision, Output = ExcitationBase>,
    SG: ElementWise<Float, Output = SGM>,
    SGM: PrecisionInput<TimePrecision, Output = SGOut>,
    SGOut: ElementWise<ExcitationBase, Output = Denominator>,  
    SGOut: ElementWise<Denominator, Output= GroundWeight>,
    ExcitationBase: ElementWise<Denominator, Output = ExcitedWeight>,

{   
    let excited_state = first_order_delocalised_rate_equation(
        e,
        s_frequency_e,
        temp,
    )?;

    let s_frequency_g_prec = s_frequency_g.element_mul(&1.0)?.map_to_precision();
    let denominator = s_frequency_g_prec.element_add(&excited_state)?;

    let ground_weight = s_frequency_g_prec.element_div(&denominator)?;
    let excited_weight = excited_state.element_div(&denominator)?;
   
    Some((ground_weight, excited_weight))
}

/// Calculate the net localised ground-state population change.
///
/// Relaxation from the excited state contributes positively, while thermal
/// excitation out of the ground state contributes negatively.
pub fn thermal_ground_state_rate
< E, SE, SG, NE, NG, Temp, ENeg, KT, Ratio, Exp, ExcitationRaw,
  ExcitationBase, PopulationRaw, PopulationOut, ExcitationOut, RelaxationRaw,
  RelaxationOut, V, >(
    e: &E,
    s_frequency_e: &SE,
    s_frequency_g: &SG,
    n_e: &NE,
    n_g: &NG,
    temp: &Temp,
) -> Option<V>
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<SE, Output = ExcitationRaw>,
    ExcitationRaw: PrecisionInput<TimePrecision, Output = ExcitationBase>,
    NG: ElementWise<Float, Output = PopulationRaw>,
    PopulationRaw: PrecisionInput<TimePrecision, Output = PopulationOut>,
    ExcitationBase: ElementWise<PopulationOut, Output = ExcitationOut>,
    SG: ElementWise<NE, Output = RelaxationRaw>,
    RelaxationRaw: PrecisionInput<TimePrecision, Output = RelaxationOut>,
    RelaxationOut: ElementWise<ExcitationOut, Output = V>,
{
    let per_particle_leaving = first_order_delocalised_rate_equation(
        e,
        s_frequency_e,
        temp,
    )?;
    let population = n_g.element_mul(&1.0)?.map_to_precision();
    let leaving = per_particle_leaving.element_mul(&population)?;
    let coming = s_frequency_g
        .element_mul(n_e)?
        .map_to_precision();

    coming.element_sub(&leaving)
}

/// Calculate the net localised excited-state population change.
///
/// Thermal excitation from the ground state and the relaxation term are
/// combined using the sign convention required by the existing simulation.
pub fn thermal_excited_state_rate
< E, SE, SG, NE, NG, Temp, ENeg, KT, Ratio, Exp, ExcitationRaw,
  ExcitationBase, PopulationRaw, PopulationOut, ExcitationOut, RelaxationRaw,
  RelaxationOut, V, >(
    e: &E,
    s_frequency_e: &SE,
    s_frequency_g: &SG,
    n_e: &NE,
    n_g: &NG,
    temp: &Temp,
) -> Option<V>
where
    E: ElementWise<Float, Output = ENeg>,
    Temp: ElementWise<Float, Output = KT>,
    ENeg: ElementWise<KT, Output = Ratio>,
    Ratio: ElementWiseUnary<Output = Exp>,
    Exp: ElementWise<SE, Output = ExcitationRaw>,
    ExcitationRaw: PrecisionInput<TimePrecision, Output = ExcitationBase>,
    NG: ElementWise<Float, Output = PopulationRaw>,
    PopulationRaw: PrecisionInput<TimePrecision, Output = PopulationOut>,
    ExcitationBase: ElementWise<PopulationOut, Output = ExcitationOut>,
    SG: ElementWise<NE, Output = RelaxationRaw>,
    RelaxationRaw: PrecisionInput<TimePrecision, Output = RelaxationOut>,
    RelaxationOut: ElementWise<ExcitationOut, Output = V>,
{
    let per_particle_coming = first_order_delocalised_rate_equation(
        e,
        s_frequency_e,
        temp,
    )?;
    let population = n_g.element_mul(&1.0)?.map_to_precision();
    let coming = per_particle_coming.element_mul(&population)?;
    let leaving = s_frequency_g
        .element_mul(n_e)?
        .map_to_precision();

    leaving.element_add(&coming)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    const TOLERANCE: f64 = 1.0e-10;

    fn expected_first_order_rate(
        e_cb: Float,
        s_frequency: Float,
        temp: Float,
    ) -> f64 {
        let exponent = (-e_cb / (BOLTZMANN_EV * temp)).exp();
        (s_frequency * exponent) as f64
    }

    fn assert_close(left: f64, right: f64) {
        assert!(
            (left - right).abs() < TOLERANCE,
            "left = {left}, right = {right}, diff = {}",
            (left - right).abs()
        );
    }

    #[test]
    fn first_order_delocalised_rate_equation_returns_expected_float_value() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 1.0e12;
        let temp: Float = 450.0;

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .expect("rate equation should return a value");

        let expected = expected_first_order_rate(e_cb, s_frequency, temp);
        assert_close(actual, expected);

        assert_eq!(e_cb, 0.45);
        assert_eq!(s_frequency, 1.0e12);
        assert_eq!(temp, 450.0);
    }

    #[test]
    fn first_order_delocalised_rate_equation_returns_expected_vec_values() {
        let e_cb: Vec<Float> = vec![0.35, 0.45, 0.55];
        let s_frequency: Vec<Float> = vec![1.0e12, 2.0e12, 3.0e12];
        let temp: Vec<Float> = vec![300.0, 450.0, 600.0];

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .expect("rate equation should return a value");

        let expected: Vec<f64> = e_cb
            .iter()
            .zip(s_frequency.iter())
            .zip(temp.iter())
            .map(|((e_cb, s_frequency), temp)| {
                expected_first_order_rate(*e_cb, *s_frequency, *temp)
            })
            .collect();

        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_close(*actual, *expected);
        }

        assert_eq!(e_cb, vec![0.35, 0.45, 0.55]);
        assert_eq!(s_frequency, vec![1.0e12, 2.0e12, 3.0e12]);
        assert_eq!(temp, vec![300.0, 450.0, 600.0]);
    }

    #[test]
    fn first_order_delocalised_rate_equation_returns_expected_ndarray_values() {
        let e_cb = array![0.35 as Float, 0.45 as Float, 0.55 as Float];
        let s_frequency = array![1.0e12 as Float, 2.0e12 as Float, 3.0e12 as Float];
        let temp = array![300.0 as Float, 450.0 as Float, 600.0 as Float];

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .expect("rate equation should return a value");

        let expected = array![
            expected_first_order_rate(e_cb[0], s_frequency[0], temp[0]),
            expected_first_order_rate(e_cb[1], s_frequency[1], temp[1]),
            expected_first_order_rate(e_cb[2], s_frequency[2], temp[2]),
        ];

        assert_eq!(actual.shape(), expected.shape());
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_close(*actual, *expected);
        }

        assert_eq!(e_cb, array![0.35 as Float, 0.45 as Float, 0.55 as Float]);
        assert_eq!(s_frequency, array![1.0e12 as Float, 2.0e12 as Float, 3.0e12 as Float]);
        assert_eq!(temp, array![300.0 as Float, 450.0 as Float, 600.0 as Float]);
    }

    fn assert_vec_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_close(*actual, *expected);
        }
    }

    #[test]
    fn exponential_energy_over_kb_t_returns_expected_values_for_all_supported_inputs() {
        let e: Float = 0.45;
        let temp: Float = 450.0;
        let scalar = exponential_energy_over_kb_t(&e, &temp)
            .expect("exponential should be calculated");
        let expected_scalar = (-e / (BOLTZMANN_EV * temp)).exp() as f64;
        assert_close(scalar as f64, expected_scalar);

        let energies: Vec<Float> = vec![0.35, 0.45, 0.55];
        let temperatures: Vec<Float> = vec![300.0, 450.0, 600.0];
        let vector = exponential_energy_over_kb_t(&energies, &temperatures)
            .expect("vector exponential should be calculated");
        let expected_vector: Vec<Float> = energies
            .iter()
            .zip(temperatures.iter())
            .map(|(energy, temperature)| {
                (-energy / (BOLTZMANN_EV * temperature)).exp()
            })
            .collect();
        assert_eq!(vector.len(), expected_vector.len());
        for (actual, expected) in vector.iter().zip(expected_vector.iter()) {
            assert_close(*actual as f64, *expected as f64);
        }

        let energies = array![0.35 as Float, 0.45 as Float, 0.55 as Float];
        let temperatures = array![300.0 as Float, 450.0 as Float, 600.0 as Float];
        let array_result = exponential_energy_over_kb_t(&energies, &temperatures)
            .expect("array exponential should be calculated");
        assert_eq!(array_result.shape(), energies.shape());
        for ((actual, energy), temperature) in array_result
            .iter()
            .zip(energies.iter())
            .zip(temperatures.iter())
        {
            let expected = (-energy / (BOLTZMANN_EV * temperature)).exp();
            assert_close(*actual as f64, expected as f64);
        }
    }

    #[test]
    fn second_order_delocalised_rate_equation_returns_per_particle_rate() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 2.0e12;
        let temp: Float = 450.0;

        let actual = second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .expect("second-order rate should be calculated");

        let expected = (s_frequency
            * (-e_cb / (BOLTZMANN_EV * temp)).exp()) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn general_order_delocalised_rate_equation_returns_per_particle_rate() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 2.0e12;
        let temp: Float = 450.0;
        let order: Float = 1.5;

        let actual = general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
            &order,
        )
        .expect("general-order rate should be calculated");

        let expected = (s_frequency
            * (-e_cb / (BOLTZMANN_EV * temp)).exp()) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn general_order_matches_first_and_second_order_at_integer_orders() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 2.0e12;
        let temp: Float = 450.0;

        let first = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .unwrap();
        let general_first = general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
            &1.0,
        )
        .unwrap();
        assert_close(first, general_first);

        let second = second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .unwrap();
        let general_second = general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
            &2.0,
        )
        .unwrap();
        assert_close(second, general_second);
    }

    #[test]
    fn hole_change_delocalised_rate_equation_returns_recombination_loss() {
        let m: Float = 0.4;
        let nc: Float = 0.3;
        let recomb: Float = 2.5;

        let actual = hole_change_delocalised_rate_equation(&m, &nc, &recomb)
            .expect("hole change should be calculated");
        let expected = (m * nc * recomb) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn retrapping_change_delocalised_rate_equation_uses_available_traps() {
        let n_tot: Float = 1.0;
        let n: Float = 0.35;
        let nc: Float = 0.2;
        let retrap: Float = 3.0;

        let actual = retrapping_change_delocalised_rate_equation(
            &n_tot,
            &n,
            &nc,
            &retrap,
        )
        .expect("retrapping change should be calculated");
        let expected = (nc * (n_tot - n) * retrap) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn trap_change_is_sum_of_thermal_release_and_retrapping() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 1.0e12;
        let n: Float = 0.25;
        let temp: Float = 450.0;
        let n_tot: Float = 1.0;
        let nc: Float = 0.1;
        let retrap: Float = 2.0;

        let actual = trap_change_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
            &n_tot,
            &nc,
            &retrap,
        )
        .expect("trap change should be calculated");
        let to_cb = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .unwrap();
        let back_cb = retrapping_change_delocalised_rate_equation(
            &n_tot,
            &n,
            &nc,
            &retrap,
        )
        .unwrap();
        assert_close(actual, back_cb - n as f64 * to_cb);
    }

    #[test]
    fn cb_band_change_combines_retrapping_recombination_and_release() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 1.0e12;
        let n: Float = 0.25;
        let temp: Float = 450.0;
        let n_tot: Float = 1.0;
        let nc: Float = 0.1;
        let m: Float = 0.2;
        let retrap: Float = 2.0;
        let recomb: Float = 3.0;

        let actual = cb_band_change_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
            &n_tot,
            &nc,
            &m,
            &retrap,
            &recomb,
        )
        .expect("conduction-band change should be calculated");

        let from_trap = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .unwrap();
        let to_trap = retrapping_change_delocalised_rate_equation(
            &n_tot,
            &n,
            &nc,
            &retrap,
        )
        .unwrap();
        let to_hole = hole_change_delocalised_rate_equation(&m, &nc, &recomb)
            .unwrap();
        assert_close(actual, to_trap - to_hole - n as f64 * from_trap);
    }

    #[test]
    fn quasi_equilibrium_rate_applies_recombination_branching_fraction() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 1.0e12;
        let n: Float = 0.25;
        let temp: Float = 450.0;
        let n_tot: Float = 1.0;
        let m: Float = 0.2;
        let retrap: Float = 2.0;
        let recomb: Float = 3.0;

        let actual = quasi_equ_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
            &n_tot,
            &m,
            &retrap,
            &recomb,
        )
        .expect("quasi-equilibrium rate should be calculated");

        let to_cb = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .unwrap();
        let to_trap = retrap * (n_tot - n);
        let to_hole = m * recomb;
        let expected = n as f64 * to_cb * (to_hole / (to_trap + to_hole)) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn tunnelling_rate_matches_exponential_distance_decay() {
        let alpha: Float = 2.0;
        let b: Float = 5.0e6;
        let r: Float = 0.4;

        let actual = tunnelling_rate(&alpha, &b, &r)
            .expect("tunnelling rate should be calculated");
        let expected = (b * (-alpha * r).exp()) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn filling_rate_scales_with_dose_rate_and_available_traps() {
        let d0: Float = 4.0;
        let d_dot: Float = 2.0;
        let n: Float = 0.25;
        let n_tot: Float = 1.0;

        let actual = filling_rate(&d0, &d_dot, &n, &n_tot)
            .expect("filling rate should be calculated");
        let expected = (d_dot / d0 * (n_tot - n)) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn filling_with_recombination_applies_retrapping_branching_fraction() {
        let d0: Float = 4.0;
        let d_dot: Float = 2.0;
        let n: Float = 0.25;
        let n_tot: Float = 1.0;
        let m: Float = 0.2;
        let retrap: Float = 2.0;
        let recomb: Float = 3.0;

        let actual = filling_with_recombination_rate(
            &d0,
            &d_dot,
            &n,
            &n_tot,
            &m,
            &retrap,
            &recomb,
        )
        .expect("filling with recombination should be calculated");

        let fill = ((d_dot / d0) * (n_tot - n) )as f64;
        let to_trap = retrap * (n_tot - n);
        let to_hole = m * recomb;
        let numer = to_trap+to_hole;
        let frac = (to_trap / numer ) as f64 ;
        let expected = fill * frac  ;
     
        assert_close(actual, expected);
    }

    #[test]
    fn affinity_probabilities_are_complementary() {
        let n: Float = 0.25;
        let n_tot: Float = 1.0;
        let m: Float = 0.5;
        let retrap: Float = 2.0;
        let recomb: Float = 1.0;

        let retrapping = retrapping_by_affinity(&n, &n_tot, &m, &retrap, &recomb)
            .expect("retrapping probability should calculate");
        let recombination = recombination_by_affinity(&n, &n_tot, &m, &retrap, &recomb)
            .expect("recombination probability should calculate");

        assert_close(retrapping, 0.75);
        assert_close(recombination, 0.25);
        assert_close(retrapping + recombination, 1.0);
    }

    #[test]
    fn affinity_probabilities_support_vector_populations() {
        let n: Vec<Float> = vec![0.0, 0.5];
        let n_tot: Float = 1.0;
        let m: Float = 0.5;
        let retrap: Float = 2.0;
        let recomb: Float = 1.0;

        let retrapping = retrapping_by_affinity(&n, &n_tot, &m, &retrap, &recomb)
            .expect("mixed scalar and vector inputs should calculate");
        let recombination = recombination_by_affinity(&n, &n_tot, &m, &retrap, &recomb)
            .expect("mixed scalar and vector inputs should calculate");

        assert_vec_close(&retrapping, &[0.8, 2.0 / 3.0]);
        assert_vec_close(&recombination, &[0.2, 1.0 / 3.0]);
    }

    #[test]
    fn distance_retrapping_factor_matches_documented_equation() {
        let prefactor: Float = 2.0;
        let mu: Float = 2.0;
        let distance: Float = 4.0;

        let actual = retrapping_probability_by_r(&prefactor, &mu, &distance)
            .expect("distance retrapping factor should calculate");
        let expected = prefactor * ((distance / mu).powi(2)).exp();

        assert_close(actual, expected);
    }

    #[test]
    fn ground_excited_state_weights_follow_boltzmann_distribution() {
        let e: Float = 0.045;
        let s_frequency_e: Float = 1.0e11;
        let s_frequency_g: Float = 1.0e9;
        let temp: Float = 450.0;

        let (ground, excited) = ground_excited_state_weights(&e, &s_frequency_e, &s_frequency_g, &temp)
            .expect("thermal state weights should calculate");
        
        let boltzmann = (-e / (BOLTZMANN_EV * temp)).exp() as f64;
        let excitation = boltzmann * s_frequency_e as f64; 
        let denominator = excitation + s_frequency_g as f64;


        let expected_ground =   s_frequency_g as f64/denominator;
        let expected_excited = excitation/denominator;
        
        assert_close(excited as f64, expected_excited as f64);
        assert_close(ground as f64, expected_ground as f64);
        assert!((ground + excited - 1.0).abs() <= (2.0 * Float::EPSILON) as f64);
    }

    #[test]
    fn ground_excited_state_weights_support_vector_energy_and_scalar_temperature() {
        let e: Vec<Float> = vec![0.01, 0.03, 0.05];
        let s_frequency_e: Vec<Float> =  vec![1.0e12,1.0e12,1.0e11];
        let s_frequency_g: Vec<Float> =  vec![1.0e9, 1.0e12,1.0e8];
        let temp: Float = 450.0;

        let (ground, excited) = ground_excited_state_weights(&e, &s_frequency_e, &s_frequency_g, &temp)
            .expect("vector thermal state weights should calculate");

        assert_eq!(ground.len(), e.len());
        assert_eq!(excited.len(), e.len());
        for (ground_weight, excited_weight) in ground.iter().zip(excited.iter()) {
            assert!((ground_weight + excited_weight - 1.0).abs() <= (2.0 * Float::EPSILON) as f64);
        }
    }

    #[test]
    fn thermal_ground_state_rate_subtracts_thermal_leaving_rate_from_arrival() {
        let e: Float = 0.45;
        let s_frequency_e: Float = 1.0e12;
        let s_frequency_g: Float = 2.0;
        let n_e: Float = 0.3;
        let n_g: Float = 0.25;
        let temp: Float = 450.0;

        let actual = thermal_ground_state_rate(
            &e,
            &s_frequency_e,
            &s_frequency_g,
            &n_e,
            &n_g,
            &temp,
        )
        .expect("ground-state rate should be calculated");

        let leaving = first_order_delocalised_rate_equation(
            &e,
            &s_frequency_e,
            &temp,
        )
        .unwrap();
        let coming = (s_frequency_g * n_e) as f64;
        assert_close(actual,  coming - n_g as f64 * leaving);
    }

    #[test]
    fn thermal_excited_state_rate_adds_thermal_arrival_to_leaving_rate() {
        let e: Float = 0.45;
        let s_frequency_e: Float = 1.0e12;
        let s_frequency_g: Float = 2.0;
        let n_e: Float = 0.3;
        let n_g: Float = 0.25;
        let temp: Float = 450.0;

        let actual = thermal_excited_state_rate(
            &e,
            &s_frequency_e,
            &s_frequency_g,
            &n_e,
            &n_g,
            &temp,
        )
        .expect("excited-state rate should be calculated");

        let coming = first_order_delocalised_rate_equation(
            &e,
            &s_frequency_e,
            &temp,
        )
        .unwrap();
        let leaving = (s_frequency_g * n_e) as f64;
        assert_close(actual, n_g as f64 * coming + leaving);
    }

    #[test]
    fn delocalised_rates_accept_vector_properties_and_scalar_temperature() {
        let e_cb: Vec<Float> = vec![0.35, 0.45, 0.55];
        let s_frequency: Vec<Float> = vec![1.0e12, 2.0e12, 3.0e12];
        let temp: Float = 450.0;

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .expect("mixed vector and scalar inputs should be supported");

        let expected: Vec<f64> = e_cb
            .iter()
            .zip(s_frequency.iter())
            .map(|(energy, frequency)| {
                (frequency * (-energy / (BOLTZMANN_EV * temp)).exp()) as f64
            })
            .collect();

        assert_vec_close(&actual, &expected);
    }

    #[test]
    fn tunnelling_rate_accepts_scalar_parameters_and_vector_distances() {
        let alpha: Float = 2.0;
        let b: Float = 5.0e6;
        let r: Vec<Float> = vec![0.1, 0.2, 0.4];

        let actual = tunnelling_rate(&alpha, &b, &r)
            .expect("mixed tunnelling inputs should be supported");

        let expected: Vec<f64> = r
            .iter()
            .map(|distance| (b * (-alpha * distance).exp()) as f64)
            .collect();

        assert_vec_close(&actual, &expected);
    }

    #[test]
    fn filling_rate_accepts_scalar_dose_and_vector_occupancy() {
        let d0: Float = 4.0;
        let d_dot: Float = 2.0;
        let n: Vec<Float> = vec![0.1, 0.3, 0.8];
        let n_tot: Float = 1.0;

        let actual = filling_rate(&d0, &d_dot, &n, &n_tot)
            .expect("mixed filling inputs should be supported");

        let expected: Vec<f64> = n
            .iter()
            .map(|occupied| (d_dot / d0 * (n_tot - occupied)) as f64)
            .collect();

        assert_vec_close(&actual, &expected);
    }

    #[test]
    fn thermal_excited_state_rate_accepts_vector_populations_and_scalar_temperature() {
        let e: Vec<Float> = vec![0.35, 0.45];
        let s_frequency_e: Vec<Float> = vec![1.0e12, 2.0e12];
        let s_frequency_g: Float = 2.0;
        let n_e: Vec<Float> = vec![0.3, 0.4];
        let n_g: Vec<Float> = vec![0.25, 0.35];
        let temp: Float = 450.0;

        let actual = thermal_excited_state_rate(
            &e,
            &s_frequency_e,
            &s_frequency_g,
            &n_e,
            &n_g,
            &temp,
        )
        .expect("mixed thermal inputs should be supported");

        assert_eq!(actual.len(), e.len());
    }
}
