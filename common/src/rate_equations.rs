/// This module contains all of the rate equations that can be used 
/// for electron transitions. They all take the following general form
/// pub fn rate_equation_name<P, A, B, C, D, E, F, G, H, I, J, K, L, T, V, W>(a, b, c) -> Option< >
/// where 
///     P: ElementWise<Float, Output = A> + ElementWise<P, Output = D> +
///        ElementWise<D, Output= H> + ElementWiseUnary<Output = E> +  
///        ElementWise<E, Output = F> + ElementWise<L, Output = V>,
///    A: ElementWise<A, Output = B> + ElementWise<P, Output = G>,   
///    B: ElementWiseUnary<Output = C>,
///    C: ElementWise<D, Output = V> + ElementWise<F, Output = V>,
///    D: ElementWise<D, Output = V>,
///    G: ElementWise<P, Output = V> + ElementWiseUnary<Output = L>,
///    H: ElementWise<P, Output = V> +  ElementWise<D, Output = I> +
///       ElementWise<Float, Output = J>,
///    I: ElementWise<J, Output = K>,
///    K: ElementWise<I, Output = W>,
///    V: PrecisionInput<TimePrecision> + PrecisionInput<TimePrecision, Output = W>,
///    W: ElementWise<W, Output = W>,
///    T: Numeric,
/// { 
///     Calculation here ...
///     Some(value.map_to_precision()?)
/// }
/// This structure allows float, Vec and ndarray to be passed to the same functions.

use crate::numeric::{Float, Numeric, ElementWise, ElementWiseUnary, PrecisionInput, TimePrecision};
use crate::constants::physical_constants::{BOLTZMANN_EV};

/// Function to calculate the exponential of the energy over KbT
pub fn exponential_energy_over_kb_t<P, A, B, C>( e: &P,temp: &P,)-> Option<C>
where
    P: ElementWise<Float, Output = A>, 
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,

{
    Some(e.element_mul(&-1.0)?
        .element_div(&temp.element_mul(&BOLTZMANN_EV)?)?
        .element_exp())
}

/// Delocalised Transitions 
/// First order rate equation 
/// e_cb conduction band energy (eV); s_frequency factor (s^-1); temp (K), n concentration of traps
pub fn first_order_delocalised_rate_equation<P, A, B, C, D, V>(
    e_cb: &P,
    s_frequency: &P,
    n: &P,
    temp: &P,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>

where
    P: ElementWise<Float, Output = A> + ElementWise<P, Output = D>, 
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<D, Output = V>,
    V: PrecisionInput<TimePrecision>,

{
    
    let exponent = exponential_energy_over_kb_t(e_cb,temp)?;

    let prefactor = s_frequency.element_mul(n)?;

    Some(exponent.element_mul(&prefactor)?
                  .multiply_to_precision(-1.0)
    )

}

/// Second order rate equation 
/// e_cb conduction band energy (eV); s_frequency factor (m^3s^-1); temp (K), n concentration of traps
pub fn second_order_delocalised_rate_equation<P, A, B, C, F, E, V>(
    e_cb: &P,
    s_frequency: &P,
    n: &P,
    temp: &P,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>

where
    P: ElementWise<Float, Output = A> + ElementWiseUnary<Output = E> + 
       ElementWise<E, Output = F>,
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<F, Output = V>,
    V: PrecisionInput<TimePrecision>,

{
    let exponent = exponential_energy_over_kb_t(e_cb,temp)?;

    let prefactor = s_frequency.element_mul(&n.element_powf(2.0))?;

    Some(exponent.element_mul(&prefactor)?
                  .multiply_to_precision(-1.0)
    )

}

/// General order rate equation 
/// e_cb conduction band energy (eV); s_frequency factor (m^3(b-1)s^-1); temp (K), n concentration of traps
pub fn general_order_delocalised_rate_equation<P, A, B, C, F, E, T, V,>(
    e_cb: &P,
    s_frequency: &P,
    n: &P,
    temp: &P,
    order: &T,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>

where
    P: ElementWise<Float, Output = A> + ElementWiseUnary<Output = E> +
       ElementWise<E, Output = F>,
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<F, Output = V>,
    T: Numeric,
    V: PrecisionInput<TimePrecision>,
{
    let exponent = exponential_energy_over_kb_t(e_cb,temp)?;

    let prefactor = s_frequency.element_mul(&n.element_powf(*order))?;

    Some(exponent.element_mul(& prefactor)?
                  .multiply_to_precision(-1.0)
    )

}

/// General order approach
/// Here the three rate equations will be given so the change in holes, traps and conduction band carriers can be calculated.
/// The quasi-equilibrium assumption leading to a single rate equation is also given

/// Change in hole concentration 
/// m concentration of holes; nc concentration of charge carriers; recomb recombination probability  
pub fn hole_change_delocalised_rate_equation<P, A, G, V>(
    m: &P,
    nc: &P,
    recomb: &P,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>

where
    P: ElementWise<Float, Output = A>,
    A: ElementWise<P, Output = G>,
    G: ElementWise<P, Output = V>,
    V: PrecisionInput<TimePrecision>,

{
    Some(
        nc.element_mul(&-1.0)?
          .element_mul(m)?
          .element_mul(recomb)?
          .map_to_precision()
    )
}

/// n concentration of traps; nc concentration of charge carriers; N_tot total number of traps; retrap recombination probability  
pub fn retrapping_change_delocalised_rate_equation<P, D, H, V>(
    n_tot: &P,
    n: &P,
    nc: &P,
    retrap: &P,
) -> Option<<V as PrecisionInput<TimePrecision>>::Output>

where
    P: ElementWise<P, Output = D> + ElementWise<D, Output= H>,
    H: ElementWise<P, Output = V>,
    V: PrecisionInput<TimePrecision>,

{   
    Some(
        nc.element_mul(&n_tot.element_sub(n)?)?
          .element_mul(retrap)?
          .map_to_precision()
    )
}
/// Change in trap concentration
/// e_cb conduction band energy (eV); s_frequency factor (s^-1); temp (K), n concentration of traps
/// n concentration of traps; nc concentration of charge carriers; N_tot total number of traps; retrap recombination probability
pub fn trap_change_delocalised_rate_equation<P, A, B, C, D, H, V, W>(
    e_cb: &P,
    s_frequency: &P,
    n: &P,
    temp: &P,
    n_tot: &P,
    nc: &P,
    retrap: &P,
) -> Option<W> 
where 
    P: ElementWise<Float, Output = A> + ElementWise<P, Output = D> +
       ElementWise<D, Output= H>,
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<D, Output = V>, 
    H: ElementWise<P, Output = V>,
    V: PrecisionInput<TimePrecision, Output = W>,
    W: ElementWise<W, Output = W>,
    
{
    let to_cb = first_order_delocalised_rate_equation(e_cb,s_frequency,n,temp)?;
    let back_cb= retrapping_change_delocalised_rate_equation(n_tot,n,nc,retrap)?;
    
    Some(to_cb.element_add(&back_cb)?)
}

/// Change in charge carrier concentration
/// e_cb conduction band energy (eV); s_frequency factor (s^-1); temp (K), n concentration of traps
/// n concentration of traps; nc concentration of charge carriers; N_tot total number of traps; retrap recombination probability
/// m concentration of holes; nc concentration of charge carriers; recomb recombination probability  
pub fn cb_band_change_delocalised_rate_equation<P, A, B, C, D, G, H, V, W>(
    e_cb: &P,
    s_frequency: &P,
    n: &P,
    temp: &P,
    n_tot: &P,
    nc: &P,
    m: &P,
    retrap: &P,
    recomb: &P,
    
) -> Option<W>

where
    P: ElementWise<Float, Output = A> + ElementWise<P, Output = D> +
       ElementWise<D, Output= H> + ElementWise<Float, Output = A>,
    A: ElementWise<A, Output = B> + ElementWise<P, Output = G>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<D, Output = V>,
    G: ElementWise<P, Output = V>,
    H: ElementWise<P, Output = V>,  
    V: PrecisionInput<TimePrecision, Output = W>,
    W: ElementWise<W, Output = W>,

{   
    let from_trap = first_order_delocalised_rate_equation(e_cb,s_frequency,n,temp)?;
    let to_trap= retrapping_change_delocalised_rate_equation(n_tot,n,nc,retrap)?;
    let to_hole = hole_change_delocalised_rate_equation(m,nc,recomb)?;
    
    Some(to_trap.element_add(&to_hole.element_sub(&from_trap)?)?)

}

/// Change in trap concentration
/// e_cb conduction band energy (eV); s_frequency factor (s^-1); temp (K), n concentration of traps
/// n concentration of traps; N_tot total number of traps; retrap recombination probability
/// m concentration of holes; nc concentration of charge carriers; recomb recombination probability  
pub fn quasi_equ_delocalised_rate_equation<P, A, B, C, D, G, H, I, J, K, V, W>(
    e_cb: &P,
    s_frequency: &P,
    n: &P,
    temp: &P,
    n_tot: &P,
    m: &P,
    retrap: &P,
    recomb: &P,
) -> Option<W> 
where 
    P: ElementWise<Float, Output = A> + ElementWise<P, Output = D> +
       ElementWise<D, Output= H>,
    A: ElementWise<A, Output = B> + ElementWise<P, Output = G>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<D, Output = V>,
    H: ElementWise<P, Output = V> +  ElementWise<D, Output = I> +
       ElementWise<Float, Output = J>,
    D: ElementWise<I, Output=K>,
    I: ElementWise<H, Output = K>,
    V: PrecisionInput<TimePrecision, Output = W>,
    W: ElementWise<K, Output = W>,
{
    let to_cb = first_order_delocalised_rate_equation(e_cb,s_frequency,n,temp)?;

    let to_trap =  retrap.element_mul(&n_tot
                            .element_sub(n)?)?;
    
    let to_hole = m.element_mul(recomb)?;
    
    let numer = to_trap.element_add(&to_hole)?;
    
    let frac = to_hole.element_div(&numer)?;
    
    Some(to_cb.element_mul(&frac)?)
}

/// Localised Transitions 
/// Function to calculate tunnelling from one state to either a hole or a trap
/// alpha exponential constant (m^-1); b attemp to escape frequency (s^-1); r distance from trap to hole (m);
pub fn tunnelling_rate<P, A, G, L, V>(
    alpha: &P,
    b: &P, 
    r: &P
)-> Option<<V as PrecisionInput<TimePrecision>>::Output>
where
    P: ElementWise<Float, Output = A> + ElementWise<L, Output = V>, 
    A: ElementWise<P, Output = G>,
    G: ElementWiseUnary<Output = L>,
    V: PrecisionInput<TimePrecision>,

{   
    Some(
        b.element_mul(&alpha
         .element_mul(&-1.0)?
         .element_mul(r)?
         .element_exp())?
         .map_to_precision()
    )
    
}

/// Filling 
pub fn filling_rate<P, D, V>(
    d0: &P,
    d_dot: &P, 
    n: &P, 
    n_tot: &P) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
 where 
    P: ElementWise<P, Output = D>,
    D: ElementWise<D, Output = V>,
    V: PrecisionInput<TimePrecision>,

 { 
    Some(
        d_dot.element_div(d0)?
             .element_mul(&n_tot
             .element_sub(n)?)?
             .map_to_precision()
    )
 }

pub fn filling_with_recombination_rate<P, D, H, I, V, W>(
    d0: &P,
    d_dot: &P,
    n: &P,
    n_tot: &P,
    m: &P,
    retrap: &P,
    recomb: &P,
) -> Option<W>
where
    P: ElementWise<P, Output = D> + ElementWise<D, Output = H>,
    D: ElementWise<D, Output = V>,
    H: ElementWise<D, Output = I> + ElementWise<I, Output = I>,
    V: PrecisionInput<TimePrecision, Output = W>,
    W: ElementWise<I, Output = W>,
{
    let fill = filling_rate(d0, d_dot, n, n_tot)?;
    
    let to_trap = retrap.element_mul(&n_tot.element_sub(n)?)?;
    let to_hole = m.element_mul(recomb)?;
    let numer = to_trap.element_add(&to_hole)?;
    let frac = to_trap.element_div(&numer)?;

    Some(fill.element_mul(&frac)?)
}

/// Internal transitions
pub fn thermal_ground_state_rate<P, A, B, C, D, V, W>(
    e: &P,
    s_frequency_e: &P,
    s_frequency_g: &P,
    n_e: &P,
    n_g: &P,
    temp: &P,) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where    
    P: ElementWise<Float, Output = A> + ElementWise<P, Output = D>, 
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<D, Output = V>,
    V: PrecisionInput<TimePrecision, Output = W>,
    W: ElementWise<D, Output = W>,
    V: PrecisionInput<TimePrecision>,
{
    let leaving  = first_order_delocalised_rate_equation(e,s_frequency_e,n_g,temp)?;
    let coming = s_frequency_g.element_mul(n_e)?;
    Some(leaving.element_add(&coming)?)

}

pub fn thermal_excited_state_rate<P, A, B, C, D, V, W>(
    e: &P,
    s_frequency_e: &P,
    s_frequency_g: &P,
    n_e: &P,
    n_g: &P,
    temp: &P,) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
where    
    P: ElementWise<Float, Output = A> + ElementWise<P, Output = D>, 
    A: ElementWise<A, Output = B>,
    B: ElementWiseUnary<Output = C>,
    C: ElementWise<D, Output = V>,
    V: PrecisionInput<TimePrecision, Output = W>,
    D: ElementWise<W, Output = W>,
    V: PrecisionInput<TimePrecision>,
{
    let coming  = first_order_delocalised_rate_equation(e,s_frequency_e,n_g,temp)?;
    let leaving = s_frequency_g.element_mul(n_e)?;
    Some(leaving.element_sub(&coming)?)

}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    const TOLERANCE: f64 = 1.0e-10;

    fn expected_first_order_rate(
        e_cb: Float,
        s_frequency: Float,
        n: Float,
        temp: Float,
    ) -> f64 {
        let exponent = (-e_cb / (BOLTZMANN_EV * temp)).exp();
        (s_frequency * n * exponent) as f64
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
        let n: Float = 0.25;
        let temp: Float = 450.0;

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .expect("rate equation should return a value");

        let expected = -1.0*expected_first_order_rate(e_cb, s_frequency, n, temp);
        assert_close(actual, expected);

        assert_eq!(e_cb, 0.45);
        assert_eq!(s_frequency, 1.0e12);
        assert_eq!(n, 0.25);
        assert_eq!(temp, 450.0);
    }

    #[test]
    fn first_order_delocalised_rate_equation_returns_expected_vec_values() {
        let e_cb: Vec<Float> = vec![0.35, 0.45, 0.55];
        let s_frequency: Vec<Float> = vec![1.0e12, 2.0e12, 3.0e12];
        let n: Vec<Float> = vec![0.1, 0.2, 0.3];
        let temp: Vec<Float> = vec![300.0, 450.0, 600.0];

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .expect("rate equation should return a value");

        let expected: Vec<f64> = e_cb
            .iter()
            .zip(s_frequency.iter())
            .zip(n.iter())
            .zip(temp.iter())
            .map(|(((e_cb, s_frequency), n), temp)| {
                -1.0*expected_first_order_rate(*e_cb, *s_frequency, *n, *temp)
            })
            .collect();

        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_close(*actual, *expected);
        }

        assert_eq!(e_cb, vec![0.35, 0.45, 0.55]);
        assert_eq!(s_frequency, vec![1.0e12, 2.0e12, 3.0e12]);
        assert_eq!(n, vec![0.1, 0.2, 0.3]);
        assert_eq!(temp, vec![300.0, 450.0, 600.0]);
    }

    #[test]
    fn first_order_delocalised_rate_equation_returns_expected_ndarray_values() {
        let e_cb = array![0.35 as Float, 0.45 as Float, 0.55 as Float];
        let s_frequency = array![1.0e12 as Float, 2.0e12 as Float, 3.0e12 as Float];
        let n = array![0.1 as Float, 0.2 as Float, 0.3 as Float];
        let temp = array![300.0 as Float, 450.0 as Float, 600.0 as Float];

        let actual = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .expect("rate equation should return a value");

        let expected = array![
            -1.0*expected_first_order_rate(e_cb[0], s_frequency[0], n[0], temp[0]),
            -1.0*expected_first_order_rate(e_cb[1], s_frequency[1], n[1], temp[1]),
            -1.0*expected_first_order_rate(e_cb[2], s_frequency[2], n[2], temp[2]),
        ];

        assert_eq!(actual.shape(), expected.shape());
        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert_close(*actual, *expected);
        }

        assert_eq!(e_cb, array![0.35 as Float, 0.45 as Float, 0.55 as Float]);
        assert_eq!(s_frequency, array![1.0e12 as Float, 2.0e12 as Float, 3.0e12 as Float]);
        assert_eq!(n, array![0.1 as Float, 0.2 as Float, 0.3 as Float]);
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
    fn second_order_delocalised_rate_equation_uses_squared_trap_concentration() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 2.0e12;
        let n: Float = 0.25;
        let temp: Float = 450.0;

        let actual = second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .expect("second-order rate should be calculated");

        let expected = (-1.0*s_frequency * n.powf(2.0)
            * (-e_cb / (BOLTZMANN_EV * temp)).exp()) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn general_order_delocalised_rate_equation_respects_the_requested_order() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 2.0e12;
        let n: Float = 0.25;
        let temp: Float = 450.0;
        let order: Float = 1.5;

        let actual = general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
            &order,
        )
        .expect("general-order rate should be calculated");

        let expected = (-1.0*s_frequency * n.powf(order)
            * (-e_cb / (BOLTZMANN_EV * temp)).exp()) as f64;
        assert_close(actual, expected);
    }

    #[test]
    fn general_order_matches_first_and_second_order_at_integer_orders() {
        let e_cb: Float = 0.45;
        let s_frequency: Float = 2.0e12;
        let n: Float = 0.25;
        let temp: Float = 450.0;

        let first = first_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .unwrap();
        let general_first = general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
            &1.0,
        )
        .unwrap();
        assert_close(first, general_first);

        let second = second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .unwrap();
        let general_second = general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
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
        let expected = -(m * nc * recomb) as f64;
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
            &n,
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
        assert_close(actual, to_cb + back_cb);
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
            &n,
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
        assert_close(actual, to_trap + to_hole - from_trap);
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
            &n,
            &temp,
        )
        .unwrap();
        let to_trap = retrap * (n_tot - n);
        let to_hole = m * recomb;
        let expected = to_cb * (to_hole / (to_trap + to_hole)) as f64;
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
    fn thermal_ground_state_rate_adds_arrival_to_thermal_leaving_rate() {
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
            &n_g,
            &temp,
        )
        .unwrap();
        let coming = (s_frequency_g * n_e) as f64;
        assert_close(actual, leaving + coming);
    }

    #[test]
    fn thermal_excited_state_rate_subtracts_thermal_arrival_from_leaving_rate() {
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
            &n_g,
            &temp,
        )
        .unwrap();
        let leaving = (s_frequency_g * n_e) as f64;
        assert_close(actual, leaving - coming);
    }

    #[test]
    fn second_order_rate_preserves_vector_and_array_shapes() {
        let e_cb: Vec<Float> = vec![0.35, 0.45, 0.55];
        let s_frequency: Vec<Float> = vec![1.0e12, 2.0e12, 3.0e12];
        let n: Vec<Float> = vec![0.1, 0.2, 0.3];
        let temp: Vec<Float> = vec![300.0, 450.0, 600.0];

        let vector = second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .unwrap();
        let expected_vector: Vec<f64> = e_cb
            .iter()
            .zip(s_frequency.iter())
            .zip(n.iter())
            .zip(temp.iter())
            .map(|(((energy, frequency), concentration), temperature)| {
                (-1.0*frequency * concentration.powf(2.0)
                    * (-energy / (BOLTZMANN_EV * temperature)).exp()) as f64
            })
            .collect();
        assert_vec_close(&vector, &expected_vector);

        let e_cb = array![0.35 as Float, 0.45 as Float, 0.55 as Float];
        let s_frequency = array![1.0e12 as Float, 2.0e12 as Float, 3.0e12 as Float];
        let n = array![0.1 as Float, 0.2 as Float, 0.3 as Float];
        let temp = array![300.0 as Float, 450.0 as Float, 600.0 as Float];
        let array_result = second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &n,
            &temp,
        )
        .unwrap();
        assert_eq!(array_result.shape(), e_cb.shape());
        for ((((actual, energy), frequency), concentration), temperature) in array_result
            .iter()
            .zip(e_cb.iter())
            .zip(s_frequency.iter())
            .zip(n.iter())
            .zip(temp.iter())
        {
            let expected = (-1.0*frequency * concentration.powf(2.0)
                * (-energy / (BOLTZMANN_EV * temperature)).exp()) as f64;
            assert_close(*actual, expected);
        }
    }
}