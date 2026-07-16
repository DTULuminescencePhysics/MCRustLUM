/// Holds the structs that are used as inputs for the various rate equations.

use crate::numeric::{ElementWise, ElementWiseUnary, Float};
use crate::rate_equations::ground_excited_state_weights;
/// Inputs shared by the ground- and excited-state delocalised equations.
///
/// `E`, `S`, and `W` independently represent the energy, frequency, and
/// weight containers. They default to [`Float`] but may be vectors or arrays.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelocalisedTransitionInputs<E = Float, S = Float, W = Float> {
    pub e_cb_ground: E,
    pub frequency_ground: S,
    pub e_cb_excited: E,
    pub frequency_excited: S,
    pub temperature: Float,
    pub ground_weight: W,
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
pub struct LocalisedTransitionInputs<Alpha = Float, B = Float, W = Float, R = Float> {
    pub alpha_ground: Alpha,
    pub frequency_ground: B,
    pub alpha_excited: Alpha,
    pub frequency_excited: B,
    pub ground_weight: W,
    pub excited_weight: W,
    pub distance: R,
}

/// Inputs for filling of the available trap population.
///
/// These fields are generic as well so a batch calculation can return the
/// same scalar, vector, or array shape for all seven transition rates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FillingTransitionInputs<D0 = Float, DDot = Float, N = Float, NTot = Float> {
    pub characteristic_dose: D0,
    pub dose_rate: DDot,
    pub occupied_population: N,
    pub total_population: NTot,
}

/// Complete scalar state needed to calculate every configured transition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransitionInputs<
    E = Float,
    S = Float,
    W = Float,
    Alpha = Float,
    B = Float,
    R = Float,
    D0 = Float,
    DDot = Float,
    N = Float,
    NTot = Float,
    Gap = Float,
> {
    pub delocalised: DelocalisedTransitionInputs<E, S, W>,
    pub localised_recombination: LocalisedTransitionInputs<Alpha, B, W, R>,
    pub localised_retrapping: LocalisedTransitionInputs<Alpha, B, W, R>,
    pub filling: FillingTransitionInputs<D0, DDot, N, NTot>,
    pub excited_energy_gap: Gap,
}

impl<E, S, W, Alpha, B, R, D0, DDot, N, NTot, Gap>
    TransitionInputs<E, S, W, Alpha, B, R, D0, DDot, N, NTot, Gap>
{
    /// Update the state-dependent inputs before calculating transition rates.
    ///
    /// The filling population is refreshed on every call. Delocalised rates
    /// are per-particle rates and therefore do not store population here. The
    /// thermal state weights are recalculated only when `temp` differs from the
    /// stored temperature, then shared by the delocalised,
    /// localised-recombination, and localised-retrapping inputs.
    pub fn update_inputs<ENeg, KT, Ratio, Exp, Denominator>(
        &mut self,
        temp: Float,
        population: N,
    ) -> Option<()>
    where
        Gap: ElementWise<Float, Output = ENeg>,
        Float: ElementWise<Float, Output = KT>
            + ElementWise<Denominator, Output = W>,
        ENeg: ElementWise<KT, Output = Ratio>,
        Ratio: ElementWiseUnary<Output = Exp>,
        Exp: ElementWise<Float, Output = Denominator>
            + ElementWise<Denominator, Output = W>,
        W: Clone,
    {
        if temp != self.delocalised.temperature {
            let (ground_weight, excited_weight) =
                ground_excited_state_weights(&self.excited_energy_gap, &temp)?;

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
        }
    }

    #[test]
    fn update_inputs_refreshes_temperature_weights_and_filling_population() {
        let mut inputs = transition_inputs();
        let new_temperature = 450.0;
        let new_population = 0.35;
        let expected_weights = ground_excited_state_weights(
            &inputs.excited_energy_gap,
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
        };
        let new_population = vec![0.4, 0.6];
        let new_temperature: Float = 450.0;
        let expected_weights = ground_excited_state_weights(
            &inputs.excited_energy_gap,
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
