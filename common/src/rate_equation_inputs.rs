//! Holds the structs that are used as inputs for the various rate equations.

use crate::numeric::{ElementWise, ElementWiseUnary, Float, TimeFloat, PrecisionInput, TimePrecision};
use crate::rate_equations::ground_excited_state_weights;

/// Identifies the physical pathway associated with a calculated rate.
///
/// Keeping the state in the identifier lets a Monte Carlo model attach a
/// distinct event to every trial lifetime. A direct solver can instead ignore
/// the identifier and sum the rates that contribute to the same population.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionKind {
    /// Ground-state tunnelling to a recombination centre.
    LocalisedRecombinationGround,
    /// Ground-state tunnelling to another trap.
    LocalisedRetrappingGround,
    /// Excited-state tunnelling to a recombination centre.
    LocalisedRecombinationExcited,
    /// Excited-state tunnelling to another trap.
    LocalisedRetrappingExcited,
    /// Ground-state release into the conduction band.
    DelocalisedGround,
    /// Excited-state release into the conduction band.
    DelocalisedExcited,
    /// Filling of an available trap.
    Filling,
}

/// One calculated rate together with the pathway that generated it.
///
/// ```
/// use common::rate_equation_inputs::{TransitionKind, TransitionRate};
///
/// let transition = TransitionRate {
///     kind: TransitionKind::Filling,
///     rate: 2.5,
/// };
/// assert_eq!(transition.rate, 2.5);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransitionRate<V = TimeFloat> {
    /// Physical event represented by this rate.
    pub kind: TransitionKind,
    /// Rate value, normally in inverse seconds.
    pub rate: V,
}

/// Inputs shared by the ground- and excited-state delocalised equations.
///
/// `E`, `S`, and `W` independently represent the energy, frequency, and
/// weight containers. They default to [`Float`] but may be vectors or arrays.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelocalisedTransitionInputs<E = Float, S = Float, W = TimeFloat> {
    /// Ground-state activation energy to the conduction band.
    pub e_cb_ground: E,
    /// Ground-state attempt-frequency factor.
    pub frequency_ground: S,
    /// Excited-state activation energy to the conduction band.
    pub e_cb_excited: E,
    /// Excited-state attempt-frequency factor.
    pub frequency_excited: S,
    /// Current temperature in kelvin.
    pub temperature: Float,
    /// Fraction of particles in the ground state.
    pub ground_weight: W,
    /// Fraction of particles in the excited state.
    pub excited_weight: W,
}

/// Inputs for one pair of ground- and excited-state tunnelling equations.
///
/// Separate instances are supplied for localised recombination and localised
/// retrapping because those pathways may use different decay constants,
/// attempt frequencies, weights, or distances.
/// `Alpha`, `B`, `W`, and `R` may each use any container combination supported
/// by the element-wise rate-equation traits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalisedTransitionInputs<Alpha = Float, B = Float, W = TimeFloat, R = Float> {
    /// Ground-state tunnelling decay constant.
    pub alpha_ground: Alpha,
    /// Ground-state tunnelling attempt frequency.
    pub frequency_ground: B,
    /// Excited-state tunnelling decay constant.
    pub alpha_excited: Alpha,
    /// Excited-state tunnelling attempt frequency.
    pub frequency_excited: B,
    /// Fraction of particles in the ground state.
    pub ground_weight: W,
    /// Fraction of particles in the excited state.
    pub excited_weight: W,
    /// Separation between the initial and target sites.
    pub distance: R,
}

/// Inputs for filling of the available trap population.
///
/// These fields are generic as well so a batch calculation can return the
/// same scalar, vector, or array shape for all seven transition rates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FillingTransitionInputs<D0 = Float, DDot = Float, N = Float, NTot = Float> {
    /// Characteristic dose controlling the filling timescale.
    pub characteristic_dose: D0,
    /// Applied dose per unit time.
    pub dose_rate: DDot,
    /// Currently occupied trap population.
    pub occupied_population: N,
    /// Total trap population available to be filled.
    pub total_population: NTot,
}

/// Complete scalar state needed to calculate every configured transition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransitionInputs<
    E = Float,
    S = Float,
    W = TimeFloat,
    Alpha = Float,
    B = Float,
    R = Float,
    D0 = Float,
    DDot = Float,
    N = Float,
    NTot = Float,
    Gap = Float,
    SE = Float,
    SG = Float,
> {
    /// Inputs for release into the conduction band.
    pub delocalised: DelocalisedTransitionInputs<E, S, W>,
    /// Inputs for localised recombination.
    pub localised_recombination: LocalisedTransitionInputs<Alpha, B, W, R>,
    /// Inputs for localised retrapping.
    pub localised_retrapping: LocalisedTransitionInputs<Alpha, B, W, R>,
    /// Inputs for trap filling.
    pub filling: FillingTransitionInputs<D0, DDot, N, NTot>,
    /// Energy gap separating the ground and excited states.
    pub excited_energy_gap: Gap,
    /// Excited-state frequency used to calculate thermal state weights.
    pub s_frequency_e: SE,
    /// Ground-state frequency used to calculate thermal state weights.
    pub s_frequency_g: SG,
}

impl<E, S, W, Alpha, B, R, D0, DDot, N, NTot, Gap, SE, SG,>
    TransitionInputs<E, S, W, Alpha, B, R, D0, DDot, N, NTot, Gap,SE, SG>
{
    /// Update the state-dependent inputs before calculating transition rates.
    ///
    /// The filling population is refreshed on every call. Delocalised rates
    /// are per-particle rates and therefore do not store population here. The
    /// thermal state weights are recalculated only when `temp` differs from the
    /// stored temperature, then shared by the delocalised,
    /// localised-recombination, and localised-retrapping inputs.
    pub fn update_inputs<ENeg, KT, Ratio, Exp, ExcitationRaw,
            ExcitationBase, Denominator, SGM, SGOut,>(
        &mut self,
        temp: Float,
        population: N,
    ) -> Option<()>
    where
        Gap: ElementWise<Float, Output = ENeg>,
        Float: ElementWise<Float, Output = KT>,
        ENeg: ElementWise<KT, Output = Ratio>,
        Ratio: ElementWiseUnary<Output = Exp>,
        Exp: ElementWise<SE, Output = ExcitationRaw>,
        ExcitationRaw: PrecisionInput<TimePrecision, Output = ExcitationBase>,
        SG: ElementWise<Float, Output = SGM>,
        SGM: PrecisionInput<TimePrecision, Output = SGOut>,
        SGOut: ElementWise<ExcitationBase, Output = Denominator>,  
        SGOut: ElementWise<Denominator, Output= W>,
        ExcitationBase: ElementWise<Denominator, Output = W>,
        W: Clone,
    {
        if temp != self.delocalised.temperature {
            let (ground_weight, excited_weight) =
                ground_excited_state_weights(&self.excited_energy_gap, &self.s_frequency_e, &self.s_frequency_g, &temp)?;

            self.delocalised.temperature = temp;
            self.delocalised.ground_weight = ground_weight.clone();
            self.delocalised.excited_weight = excited_weight.clone();

            self.localised_recombination.ground_weight = ground_weight.clone();
            self.localised_recombination.excited_weight = excited_weight.clone();

            self.localised_retrapping.ground_weight = ground_weight;
            self.localised_retrapping.excited_weight = excited_weight;
        }

        self.filling.occupied_population = population;

        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn transition_inputs() -> TransitionInputs {
        TransitionInputs {
            delocalised: DelocalisedTransitionInputs {
                e_cb_ground: 0.35,
                frequency_ground: 1.0e12,
                e_cb_excited: 0.25,
                frequency_excited: 2.0e12,
                temperature: 300.0,
                ground_weight: 0.6,
                excited_weight: 0.4,
            },
            localised_recombination: LocalisedTransitionInputs {
                alpha_ground: 1.0,
                frequency_ground: 2.0,
                alpha_excited: 1.5,
                frequency_excited: 3.0,
                ground_weight: 0.6,
                excited_weight: 0.4,
                distance: 0.5,
            },
            localised_retrapping: LocalisedTransitionInputs {
                alpha_ground: 0.8,
                frequency_ground: 4.0,
                alpha_excited: 1.2,
                frequency_excited: 5.0,
                ground_weight: 0.6,
                excited_weight: 0.4,
                distance: 0.25,
            },
            filling: FillingTransitionInputs {
                characteristic_dose: 4.0,
                dose_rate: 2.0,
                occupied_population: 0.2,
                total_population: 1.0,
            },
            excited_energy_gap: 0.045,
            s_frequency_e: 1.0e12,
            s_frequency_g: 1.0e12,
        }
    }

    #[test]
    fn update_inputs_refreshes_temperature_weights_and_filling_population() {
        let mut inputs = transition_inputs();
        let new_temperature = 450.0;
        let new_population = 0.35;
        let expected_weights = ground_excited_state_weights(
            &inputs.excited_energy_gap,
            &inputs.s_frequency_e, 
            &inputs.s_frequency_g,
            &new_temperature,
        )
        .unwrap();

        inputs
            .update_inputs(new_temperature, new_population)
            .expect("scalar transition inputs should update");

        assert_eq!(inputs.delocalised.temperature, new_temperature);
        assert_eq!(inputs.filling.occupied_population, new_population);
        assert_eq!(
            (inputs.delocalised.ground_weight, inputs.delocalised.excited_weight),
            expected_weights,
        );
        assert_eq!(
            (
                inputs.localised_recombination.ground_weight,
                inputs.localised_recombination.excited_weight,
            ),
            expected_weights,
        );
        assert_eq!(
            (
                inputs.localised_retrapping.ground_weight,
                inputs.localised_retrapping.excited_weight,
            ),
            expected_weights,
        );
    }

    #[test]
    fn update_inputs_keeps_weights_when_temperature_is_unchanged() {
        let mut inputs = transition_inputs();
        let original_delocalised_weights = (
            inputs.delocalised.ground_weight,
            inputs.delocalised.excited_weight,
        );

        inputs
            .update_inputs(inputs.delocalised.temperature, 0.75)
            .expect("population-only update should succeed");

        assert_eq!(
            (inputs.delocalised.ground_weight, inputs.delocalised.excited_weight),
            original_delocalised_weights,
        );
        assert_eq!(inputs.filling.occupied_population, 0.75);
    }

    #[test]
    fn update_inputs_supports_vector_weights_and_population() {
        let mut inputs = TransitionInputs {
            delocalised: DelocalisedTransitionInputs {
                e_cb_ground: vec![0.35, 0.45],
                frequency_ground: vec![1.0e12, 2.0e12],
                e_cb_excited: vec![0.25, 0.35],
                frequency_excited: vec![2.0e12, 3.0e12],
                temperature: 300.0,
                ground_weight: vec![0.5, 0.5],
                excited_weight: vec![0.5, 0.5],
            },
            localised_recombination: LocalisedTransitionInputs {
                alpha_ground: 1.0,
                frequency_ground: 2.0,
                alpha_excited: 1.5,
                frequency_excited: 3.0,
                ground_weight: vec![0.5, 0.5],
                excited_weight: vec![0.5, 0.5],
                distance: vec![0.25, 0.5],
            },
            localised_retrapping: LocalisedTransitionInputs {
                alpha_ground: 0.8,
                frequency_ground: 4.0,
                alpha_excited: 1.2,
                frequency_excited: 5.0,
                ground_weight: vec![0.5, 0.5],
                excited_weight: vec![0.5, 0.5],
                distance: vec![0.25, 0.5],
            },
            filling: FillingTransitionInputs {
                characteristic_dose: 4.0,
                dose_rate: 2.0,
                occupied_population: vec![0.2, 0.3],
                total_population: 1.0,
            },
            excited_energy_gap: vec![0.04, 0.05],
            s_frequency_e: vec![1.0e12, 1.0e9],
            s_frequency_g: vec![1.0e12, 1.0e9],
        };
        let new_population = vec![0.4, 0.6];
        let new_temperature: Float = 450.0;
        let expected_weights = ground_excited_state_weights(
            &inputs.excited_energy_gap,
            &inputs.s_frequency_e, 
            &inputs.s_frequency_g,
            &new_temperature,
        )
        .unwrap();

        inputs
            .update_inputs(new_temperature, new_population.clone())
            .expect("vector transition inputs should update");

        assert_eq!(inputs.delocalised.ground_weight, expected_weights.0);
        assert_eq!(inputs.delocalised.excited_weight, expected_weights.1);
        assert_eq!(inputs.filling.occupied_population, new_population);
    }
}
