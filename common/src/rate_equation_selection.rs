//! Runtime configuration for the rate equations used by the simulation.
//!
//! This module separates two choices that are often supplied by configuration:
//! which physical equation to evaluate and which transition pathways are active.
//! The selection types then dispatch to the generic implementations in
//! [`crate::rate_equations`].

use crate::numeric::{Float, ElementWise, ElementWiseUnary, PrecisionInput, TimePrecision};
use crate::rate_equations;
use crate::rate_equation_inputs::{
    DelocalisedTransitionInputs,
    FillingTransitionInputs,
    LocalisedTransitionInputs,
    TransitionInputs,
    TransitionKind,
    TransitionRate,
};
use std::str::FromStr;


/// Runtime-selectable delocalised rate equation.
/// The selected equation can be stored in a configuration structure and
/// called through the same `calculate` method. `GeneralOrder` stores the
/// kinetic order required by the general-order equation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DelocalisedRateEquationType {
    FirstOrder,
    SecondOrder,
    GeneralOrder { order: Float },
}

impl DelocalisedRateEquationType {
    /// Evaluate the selected kinetic model using a common interface.
    ///
    /// The element-wise traits allow scalar, vector, and array inputs to be
    /// mixed when the corresponding broadcasting implementation exists.
    /// `None` indicates that an element-wise operation could not be completed,
    /// for example because two input containers have incompatible shapes.
    pub fn calculate<E, S, ENeg, KT, Ratio, Exp, V>(
        &self,
        e_cb: &E,
        s_frequency: &S,
        temperature: &Float,
    ) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
    where
        E: ElementWise<Float, Output = ENeg>,
        Float: ElementWise<Float, Output = KT>,
        ENeg: ElementWise<KT, Output = Ratio>,
        Ratio: ElementWiseUnary<Output = Exp>,
        Exp: ElementWise<S, Output = V>,
        V: PrecisionInput<TimePrecision>,
    {
        match *self {
            Self::FirstOrder => rate_equations::first_order_delocalised_rate_equation(
                e_cb,
                s_frequency,
                temperature,
            ),
            Self::SecondOrder => rate_equations::second_order_delocalised_rate_equation(
                e_cb,
                s_frequency,
                temperature,
            ),
            Self::GeneralOrder { order } => {
                rate_equations::general_order_delocalised_rate_equation(
                    e_cb,
                    s_frequency,
                    temperature,
                    &order,
                )
            }
        }
    }
}

impl FromStr for DelocalisedRateEquationType {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim().to_ascii_lowercase();

        match value.as_str() {
            "none" | "off" | "disabled" => Err(
                "none selects the transition state, not a delocalised rate equation type"
                    .to_string(),
            ),
            "first" | "first_order" | "first-order" => Ok(Self::FirstOrder),
            "second" | "second_order" | "second-order" => Ok(Self::SecondOrder),
            _ => {
                if let Some(order) = value.strip_prefix("general:") {
                    let order = order.parse::<Float>().map_err(|error| {
                        format!("invalid general-order value '{order}': {error}")
                    })?;

                    if !order.is_finite() || order <= 0.0 {
                        return Err("general order must be finite and greater than zero".to_string());
                    }

                    Ok(Self::GeneralOrder { order })
                } else {
                    Err(format!(
                        "unknown delocalised rate equation '{value}'; expected first, second, or general:<order>"
                    ))
                }
            }
        }
    }
}

/// Selects the states from which delocalised release is enabled.
///
/// Active states use the embedded kinetic model. Inactive states still return
/// a zero with the same shape as that state's energy input, so callers can
/// combine both outputs without special-casing the selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DelocalisedRateEquation {
    /// Evaluate only the ground-state release rate.
    Ground{re: DelocalisedRateEquationType},
    /// Evaluate only the excited-state release rate.
    Excited{re: DelocalisedRateEquationType},
    /// Evaluate both ground- and excited-state release rates.
    Both{re: DelocalisedRateEquationType},
    /// Disable delocalised release from both states.
    None,
}

impl DelocalisedRateEquation {
    /// Calculate the weighted `(ground, excited)` release rates.
    ///
    /// Each tuple entry preserves the shape associated with its state. An
    /// inactive state produces a shape-preserving zero; `None` is reserved for
    /// failed element-wise operations such as incompatible input shapes.
    pub fn calculate<E, S, W, ENeg, KT, Ratio, Exp, WeightRaw, Weight,
        Weighted, RateRaw, V,
    >(
        &self,
        inputs: &DelocalisedTransitionInputs<E, S, W>,
    ) -> (Option<V>, Option<V>)
    where
        E: ElementWise<Float, Output = ENeg>,
        Float: ElementWise<Float, Output = KT>,
        ENeg: ElementWise<KT, Output = Ratio>
            + PrecisionInput<TimePrecision, Output = V>,
        Ratio: ElementWiseUnary<Output = Exp>,
        Exp: ElementWise<S, Output = RateRaw>,
        RateRaw: PrecisionInput<TimePrecision, Output = V>,
        W: ElementWise<Float, Output = WeightRaw>,
        WeightRaw: PrecisionInput<TimePrecision, Output = Weight>,
        Weight: ElementWise<V, Output = Weighted>,
        Weighted: PrecisionInput<TimePrecision, Output = V>,
    {
        let zero_ground = || {
            inputs.e_cb_ground
                .element_mul(&0.0)
                .map(|zero| zero.map_to_precision())
        };

        let zero_excited = || {
            inputs.e_cb_excited
                .element_mul(&0.0)
                .map(|zero| zero.map_to_precision())
        };

        let weighted_ground = |re: DelocalisedRateEquationType| {
            re.calculate(
                &inputs.e_cb_ground,
                &inputs.frequency_ground,
                &inputs.temperature,
            )
            .and_then(|rate| {
                    inputs.ground_weight.element_mul(&1.0)
                        .map(|weight| weight.map_to_precision())?
                        .element_mul(&rate)
                        .map(|weighted| weighted.map_to_precision())
                })
        };

        let weighted_excited = |re: DelocalisedRateEquationType| {
            re.calculate(
                &inputs.e_cb_excited,
                &inputs.frequency_excited,
                &inputs.temperature,
            )
            .and_then(|rate| {
                    inputs.excited_weight.element_mul(&1.0)
                        .map(|weight| weight.map_to_precision())?
                        .element_mul(&rate)
                        .map(|weighted| weighted.map_to_precision())
                })
        };

        match *self {
            Self::Ground { re } => (weighted_ground(re), zero_excited()),
            Self::Excited { re } => (zero_ground(), weighted_excited(re)),
            Self::Both { re } => (weighted_ground(re), weighted_excited(re)),
            Self::None => (zero_ground(), zero_excited()),
        }
    }
}
impl DelocalisedRateEquation {
    /// Build a selection from separate state and kinetic-model settings.
    ///
    /// `value` accepts `ground`, `excited`, `both`, or `none`; `value2`
    /// accepts `first`, `second`, or `general:<order>`. The equation setting is
    /// ignored when the state selection is `none`.
    pub fn from_strs(value: &str, value2: &str) -> Result<Self, String> {
        let selection = value.trim().to_ascii_lowercase();

        if matches!(selection.as_str(), "none" | "off" | "disabled") {
            return Ok(Self::None);
        }

        let rate_equation = value2.parse::<DelocalisedRateEquationType>()?;

        match selection.as_str() {
            "ground" | "ground_only" | "ground-only" => {
                Ok(Self::Ground { re: rate_equation })
            }
            "excited" | "excited_only" | "excited-only" => {
                Ok(Self::Excited { re: rate_equation })
            }
            "both" | "ground_excited" | "ground-excited" => {
                Ok(Self::Both { re: rate_equation })
            }
            _ => Err(format!(
                "unknown delocalised transition selection '{selection}'; expected ground, excited, both, or none"
            )),
        }
    }
}

impl FromStr for DelocalisedRateEquation {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim();
        let normalized = value.to_ascii_lowercase();

        if matches!(normalized.as_str(), "none" | "off" | "disabled") {
            return Ok(Self::None);
        }

        // Split only once so a nested value such as `general:1.5` remains
        // intact for `DelocalisedRateEquationType` to parse.
        let (selection, equation) = value.split_once(':').ok_or_else(|| {
            format!(
                "invalid delocalised selection '{value}'; expected <ground|excited|both>:<first|second|general:order> or none"
            )
        })?;

        Self::from_strs(selection, equation)
    }
}

/// Selects the states from which localised tunnelling is enabled.
///
/// As with delocalised selection, disabled states produce shape-preserving
/// zero values rather than being omitted from the returned pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LocalisedRateEquation {
    Ground,
    Excited,
    Both,
    None,
}

impl LocalisedRateEquation {
    /// Calculate the weighted `(ground, excited)` tunnelling rates.
    ///
    /// The distance input `r` defines the output shape used for inactive
    /// states. A `None` tuple entry indicates a failed element-wise operation.
    pub fn calculate
    < Alpha, B, W, R, NegativeAlpha, Product, Exp, Zero, WeightRaw,
      Weight, Weighted, RateOut, V, >(
        &self,
        inputs: &LocalisedTransitionInputs<Alpha, B, W, R>,
    ) -> (Option<V>, Option<V>)
       
    where
        Alpha: ElementWise<Float, Output = NegativeAlpha>,
        NegativeAlpha: ElementWise<R, Output = Product>,
        Product: ElementWiseUnary<Output = Exp>,
        B: ElementWise<Exp, Output = RateOut>,
        RateOut: PrecisionInput<TimePrecision,Output = V>, 
        R: ElementWise<Float, Output = Zero>,
        Zero: PrecisionInput<TimePrecision,Output = V>,
        V: PrecisionInput<TimePrecision>,
        W: ElementWise<Float, Output = WeightRaw>,
        WeightRaw: PrecisionInput<TimePrecision, Output = Weight>,
        Weight: ElementWise<V, Output = Weighted>,
        Weighted: PrecisionInput<TimePrecision, Output = V>,
    {
        let zero = || {
            inputs.distance.element_mul(&0.0)
                .map(|zero| zero.map_to_precision())
                
        };
        let weighted_ground = || {
            let rate = rate_equations::tunnelling_rate(
                &inputs.alpha_ground,
                &inputs.frequency_ground,
                &inputs.distance,
            )?;
            Some(
                inputs.ground_weight.element_mul(&1.0)?
                .map_to_precision()
                .element_mul(&rate)?
                .map_to_precision())

        };
        let weighted_excited = || {
            let rate = rate_equations::tunnelling_rate(
                &inputs.alpha_excited,
                &inputs.frequency_excited,
                &inputs.distance,
            )?;
            Some(
                inputs.excited_weight.element_mul(&1.0)?
                .map_to_precision()
                .element_mul(&rate)?
                .map_to_precision())
        };

        match *self {
            Self::Ground => (weighted_ground(), zero()),
            Self::Excited => (zero(), weighted_excited()),
            Self::Both => (weighted_ground(), weighted_excited()),
            Self::None => (zero(), zero()),
        }
    }
}
impl FromStr for LocalisedRateEquation {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim().to_ascii_lowercase();

        match value.as_str() {
            "none" | "off" | "disabled" => Ok(Self::None),
            "ground" | "ground_only" | "ground-only" => Ok(Self::Ground),
            "excited" | "excited_only" | "excited-only" => Ok(Self::Excited),
            "both" | "ground_excited" | "ground-excited" => Ok(Self::Both),
            _ => { 
                Err(format!(
                        "unknown localised rate equation '{value}'; expected ground, excited, both or none"
                    ))
            }
        }
    }
}


/// Runtime selection for filling transitions.
/// Add further filling models as new enum variants. Variants can carry
/// model-specific configuration values when needed.

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FillingRateEquation {
    Basic,
    None, 
}

impl FillingRateEquation {
    /// Evaluate filling, or return a zero matching the characteristic-dose
    /// shape when filling is disabled.
    pub fn calculate<
        D0, DDot, N, NTot, DoseRatio, Available, Zero, V, >(
        &self,
        inputs: &FillingTransitionInputs<D0, DDot, N, NTot>,
    ) -> Option<<V as PrecisionInput<TimePrecision>>::Output>
    where
        DDot: ElementWise<D0, Output = DoseRatio>,
        NTot: ElementWise<N, Output = Available>,
        DoseRatio: ElementWise<Available, Output = V>,
        D0: ElementWise<Float, Output = Zero>,
        Zero: PrecisionInput<
            TimePrecision,
            Output = <V as PrecisionInput<TimePrecision>>::Output,
        >,
        V: PrecisionInput<TimePrecision>,
    {
        match *self {
            Self::Basic => rate_equations::filling_rate(
                &inputs.characteristic_dose,
                &inputs.dose_rate,
                &inputs.occupied_population,
                &inputs.total_population,
            ),
            Self::None => Some(
                inputs.characteristic_dose.element_mul(&0.0)?
                    .map_to_precision(),
            ),
        }
    }
}

impl FromStr for FillingRateEquation {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim().to_ascii_lowercase();

        match value.as_str() {
            "basic" | "fill" | "filling" => Ok(Self::Basic),
            "none" | "off" | "disabled" => Ok(Self::None),
            _ => Err(format!(
                "unknown filling rate equation '{value}'; expected basic or none"
            )),
        }
    }
}

/// All equation selections needed to evaluate a transition configuration.
///
/// `localised_recomb` is configured directly, while `localised_retrap` is
/// derived from the ground/excited flags in [`RetrappingSelection`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransitionsTypes {
    pub delocalised: DelocalisedRateEquation,
    pub localised_recomb: LocalisedRateEquation,
    pub localised_retrap: LocalisedRateEquation,
    pub filling: FillingRateEquation,
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Groups transition configurations by the two non-localised retrapping flags.
///
/// The contained [`TransitionsTypes`] retains the complete equation selection;
/// the outer variant makes conduction-band and filling retrapping cheap to
/// dispatch on in the simulation.
pub enum Transitions{
    NoCbFillRetrapping{transitions:TransitionsTypes},
    FillRetrapping{transitions:TransitionsTypes}, 
    CbRetrapping{transitions:TransitionsTypes},
    FillCbRetrapping{transitions:TransitionsTypes}
}


#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
/// Independent switches for the available retrapping pathways.
///
/// Text configuration uses the abbreviations `cb` (conduction band), `fi`
/// (filling), `gs` (ground state), and `es` (excited state). Components may be
/// compact or separated with `-` and `_`, for example `cbfigs` and `cb_fi_gs`
/// are equivalent.
pub struct RetrappingSelection {
    pub cb: bool,
    pub filling: bool,
    pub ground: bool,
    pub excited: bool,
}

impl RetrappingSelection {
    /// Return a selection with every retrapping pathway disabled.
    pub const fn none() -> Self {
        Self {
            cb: false,
            filling: false,
            ground: false,
            excited: false,
        }
    }

    /// Collapse the two state flags into the corresponding localised mode.
    pub const fn localised_rate_equation(&self) -> LocalisedRateEquation {
        match (self.ground, self.excited) {
            (true, true) => LocalisedRateEquation::Both,
            (true, false) => LocalisedRateEquation::Ground,
            (false, true) => LocalisedRateEquation::Excited,
            (false, false) => LocalisedRateEquation::None,
        }
    }

    /// Choose the outer transition variant from the CB and filling flags.
    pub const fn transitionselection(&self, transitions: TransitionsTypes) -> Transitions {
        match (self.cb, self.filling) {
            (true, true) => Transitions::FillCbRetrapping{transitions}, 
            (true, false) => Transitions::CbRetrapping{transitions},
            (false, true) => Transitions::FillRetrapping{transitions},
            (false, false) => Transitions::NoCbFillRetrapping{transitions}
        }
    }
}

impl FromStr for RetrappingSelection {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let normalized = value
            .trim()
            .to_ascii_lowercase()
            .replace(['-', '_'], "");

        if matches!(normalized.as_str(), "none" | "off" | "disabled" | "noretrapping") {
            return Ok(Self::none());
        }

        if normalized.is_empty() {
            return Err("retrapping selection cannot be empty".to_string());
        }

        let mut selection = Self::none();
        let mut remaining = normalized.as_str();

        // Consume fixed two-character tokens so both compact and separated
        // forms share one parser after separators have been removed.
        while !remaining.is_empty() {
            if let Some(rest) = remaining.strip_prefix("cb") {
                if selection.cb {
                    return Err(format!(
                        "duplicate conduction-band retrapping component in '{value}'"
                    ));
                }
                selection.cb = true;
                remaining = rest;
            } else if let Some(rest) = remaining.strip_prefix("fi") {
                if selection.filling {
                    return Err(format!(
                        "duplicate filling retrapping component in '{value}'"
                    ));
                }
                selection.filling = true;
                remaining = rest;
            } else if let Some(rest) = remaining.strip_prefix("gs") {
                if selection.ground {
                    return Err(format!(
                        "duplicate ground-state retrapping component in '{value}'"
                    ));
                }
                selection.ground = true;
                remaining = rest;
            } else if let Some(rest) = remaining.strip_prefix("es") {
                if selection.excited {
                    return Err(format!(
                        "duplicate excited-state retrapping component in '{value}'"
                    ));
                }
                selection.excited = true;
                remaining = rest;
            } else {
                return Err(format!(
                    "unknown retrapping component in '{value}' near '{remaining}'"
                ));
            }
        }

        Ok(selection)
    }
}

impl Transitions {
    /// Calculate the seven physical transitions used by the simulation.
    ///
    /// The returned array always has the same order: ground recombination,
    /// ground retrapping, excited recombination, excited retrapping, ground
    /// delocalisation, excited delocalisation, and filling. A disabled pathway
    /// occupies its normal position with a zero rate.
    ///
    /// Conduction-band retrapping is intentionally not a separate rate. The
    /// outer `Transitions` variant records whether it should be applied later
    /// by a direct solver or after the Monte Carlo model selects an event.
    ///
    /// `None` means one of the underlying element-wise operations could not
    /// produce an output, usually because input container shapes differ.
    pub fn calculate<
        E, S, W, Alpha, B, R, D0, DDot, N, NTot, Gap, SE,SG, KT, ENeg, DRatio,
        DExp, DRaw, NegativeAlpha, LProduct, LExp, LRaw, LZero, WeightRaw,
        Weight, WeightedRaw, DoseRatio, Available, FillRaw, FillZero, V,
    >(
        &self,
        inputs: &TransitionInputs<E, S, W, Alpha, B, R, D0, DDot, N, NTot, Gap, SE, SG>,
    ) -> Option<[TransitionRate<V>; 7]>
    where
        Float: ElementWise<Float, Output = KT>,
        E: ElementWise<Float, Output = ENeg>,
        ENeg: ElementWise<KT, Output = DRatio>
            + PrecisionInput<TimePrecision, Output = V>,
        DRatio: ElementWiseUnary<Output = DExp>,
        DExp: ElementWise<S, Output = DRaw>,
        DRaw: PrecisionInput<TimePrecision, Output = V>,
        Alpha: ElementWise<Float, Output = NegativeAlpha>,
        NegativeAlpha: ElementWise<R, Output = LProduct>,
        LProduct: ElementWiseUnary<Output = LExp>,
        B: ElementWise<LExp, Output = LRaw>,
        LRaw: PrecisionInput<TimePrecision, Output = V>,
        R: ElementWise<Float, Output = LZero>,
        LZero: PrecisionInput<TimePrecision, Output = V>,
        W: ElementWise<Float, Output = WeightRaw>,
        WeightRaw: PrecisionInput<TimePrecision, Output = Weight>,
        Weight: ElementWise<V, Output = WeightedRaw>,
        WeightedRaw: PrecisionInput<TimePrecision, Output = V>,
        DDot: ElementWise<D0, Output = DoseRatio>,
        NTot: ElementWise<N, Output = Available>,
        DoseRatio: ElementWise<Available, Output = FillRaw>,
        FillRaw: PrecisionInput<TimePrecision, Output = V>,
        D0: ElementWise<Float, Output = FillZero>,
        FillZero: PrecisionInput<TimePrecision, Output = V>,
        V: PrecisionInput<TimePrecision>,
    {
        let transitions = match self {
            Self::NoCbFillRetrapping { transitions }
            | Self::FillRetrapping { transitions }
            | Self::CbRetrapping { transitions }
            | Self::FillCbRetrapping { transitions } => transitions,
        };

        let (delocalised_ground, delocalised_excited) =
            transitions.delocalised.calculate(&inputs.delocalised);

        let (recombination_ground, recombination_excited) = transitions
            .localised_recomb
            .calculate(&inputs.localised_recombination);

        let (retrapping_ground, retrapping_excited) = transitions
            .localised_retrap
            .calculate(&inputs.localised_retrapping);

        let filling_rate = transitions.filling.calculate(&inputs.filling)?;

        Some([
            TransitionRate {
                kind: TransitionKind::LocalisedRecombinationGround,
                rate: recombination_ground?,
            },
            TransitionRate {
                kind: TransitionKind::LocalisedRetrappingGround,
                rate: retrapping_ground?,
            },
            TransitionRate {
                kind: TransitionKind::LocalisedRecombinationExcited,
                rate: recombination_excited?,
            },
            TransitionRate {
                kind: TransitionKind::LocalisedRetrappingExcited,
                rate: retrapping_excited?,
            },
            TransitionRate {
                kind: TransitionKind::DelocalisedGround,
                rate: delocalised_ground?,
            },
            TransitionRate {
                kind: TransitionKind::DelocalisedExcited,
                rate: delocalised_excited?,
            },
            TransitionRate {
                kind: TransitionKind::Filling,
                rate: filling_rate,
            },
        ])
    }

    /// Parse and validate the complete transition configuration.
    ///
    /// Retrapping pathways that depend on a rate equation cannot be enabled
    /// while that equation is disabled. The localised retrapping mode is
    /// derived from the `gs` and `es` components of `retrapping`.
    pub fn from_strs(
        retrapping: &str,
        delocalised_selection: &str,
        delocalised_type: &str,
        localised_recombination: &str,
        filling: &str,
    ) -> Result<Self, String> {
        let retrapping = retrapping.parse::<RetrappingSelection>()?;
        let delocalised = DelocalisedRateEquation::from_strs(
            delocalised_selection,
            delocalised_type,
        )?;
        let localised_recomb = localised_recombination
            .parse::<LocalisedRateEquation>()?;
        let filling = filling.parse::<FillingRateEquation>()?;

        if retrapping.cb && matches!(delocalised, DelocalisedRateEquation::None) {
            return Err(
                "conduction-band retrapping requires a non-none delocalised rate equation"
                    .to_string(),
            );
        }

        if retrapping.filling && matches!(filling, FillingRateEquation::None) {
            return Err(
                "filling retrapping requires a non-none filling rate equation"
                    .to_string(),
            );
        }
        let localised_retrap = retrapping.localised_rate_equation();
        let transitions = TransitionsTypes{
            delocalised,
            localised_recomb,
            localised_retrap,
            filling
        };
        Ok(retrapping.transitionselection(transitions)) 

    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::rate_equation_inputs::{
        DelocalisedTransitionInputs,
        FillingTransitionInputs,
        LocalisedTransitionInputs,
    };
    use crate::numeric::TimeFloat;

    fn scalar_transition_inputs() -> TransitionInputs {
        TransitionInputs {
            delocalised: DelocalisedTransitionInputs {
                e_cb_ground: 0.35,
                frequency_ground: 1.0e12,
                e_cb_excited: 0.25,
                frequency_excited: 2.0e12,
                temperature: 450.0,
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
                ground_weight: 0.7,
                excited_weight: 0.3,
                distance: 0.25,
            },
            filling: FillingTransitionInputs {
                characteristic_dose: 4.0,
                dose_rate: 2.0,
                occupied_population: 0.25,
                total_population: 1.0,
            },
            excited_energy_gap: 0.045,
            s_frequency_e: 1.0e12,
            s_frequency_g: 1.0e12,
        }
    }

    #[test]
    fn calculate_returns_seven_rates_with_stable_transition_kinds() {
        let selected = Transitions::FillCbRetrapping {
            transitions: TransitionsTypes {
                delocalised: DelocalisedRateEquation::Both {
                    re: DelocalisedRateEquationType::FirstOrder,
                },
                localised_recomb: LocalisedRateEquation::Both,
                localised_retrap: LocalisedRateEquation::Both,
                filling: FillingRateEquation::Basic,
            },
        };

        let rates = selected
            .calculate(&scalar_transition_inputs())
            .expect("all scalar transition rates should calculate");
        let kinds: Vec<_> = rates.iter().map(|output| output.kind).collect();

        assert_eq!(
            kinds,
            vec![
                TransitionKind::LocalisedRecombinationGround,
                TransitionKind::LocalisedRetrappingGround,
                TransitionKind::LocalisedRecombinationExcited,
                TransitionKind::LocalisedRetrappingExcited,
                TransitionKind::DelocalisedGround,
                TransitionKind::DelocalisedExcited,
                TransitionKind::Filling,
            ],
        );
        assert!(rates.iter().all(|output| output.rate.is_finite()));
    }

    #[test]
    fn calculate_represents_disabled_pathways_with_zero_rates() {
        let selected = Transitions::NoCbFillRetrapping {
            transitions: TransitionsTypes {
                delocalised: DelocalisedRateEquation::Ground {
                    re: DelocalisedRateEquationType::FirstOrder,
                },
                localised_recomb: LocalisedRateEquation::Excited,
                localised_retrap: LocalisedRateEquation::None,
                // Filling is selected independently of later CB retrapping.
                filling: FillingRateEquation::Basic,
            },
        };

        let rates = selected
            .calculate(&scalar_transition_inputs())
            .expect("selected scalar transition rates should calculate");
        assert_eq!(rates.len(), 7);
        assert_eq!(rates[0].rate, 0.0); // Ground recombination is disabled.
        assert_eq!(rates[1].rate, 0.0); // Ground retrapping is disabled.
        assert!(rates[2].rate > 0.0); // Excited recombination is enabled.
        assert_eq!(rates[3].rate, 0.0); // Excited retrapping is disabled.
        assert!(rates[4].rate > 0.0); // Ground delocalisation is enabled.
        assert_eq!(rates[5].rate, 0.0); // Excited delocalisation is disabled.
        assert!(rates[6].rate > 0.0); // Filling is controlled by its own selection.
    }

    #[test]
    fn calculate_supports_vector_transition_inputs() {
        let selected = Transitions::FillCbRetrapping {
            transitions: TransitionsTypes {
                delocalised: DelocalisedRateEquation::Both {
                    re: DelocalisedRateEquationType::FirstOrder,
                },
                localised_recomb: LocalisedRateEquation::Both,
                localised_retrap: LocalisedRateEquation::Both,
                filling: FillingRateEquation::Basic,
            },
        };
        let mut inputs = TransitionInputs {
            delocalised: DelocalisedTransitionInputs {
                e_cb_ground: vec![0.35 as Float, 0.45],
                frequency_ground: vec![1.0e12 as Float, 2.0e12],
                e_cb_excited: vec![0.25 as Float, 0.35],
                frequency_excited: vec![2.0e12 as Float, 3.0e12],
                temperature: 300.0,
                ground_weight: vec![0.5 as TimeFloat, 0.5],
                excited_weight: vec![0.5 as TimeFloat, 0.5],
            },
            localised_recombination: LocalisedTransitionInputs {
                alpha_ground: 1.0 as Float,
                frequency_ground: 2.0 as Float,
                alpha_excited: 1.5 as Float,
                frequency_excited: 3.0 as Float,
                ground_weight: vec![0.5 as TimeFloat, 0.5],
                excited_weight: vec![0.5 as TimeFloat, 0.5],
                distance: vec![0.25 as Float, 0.5],
            },
            localised_retrapping: LocalisedTransitionInputs {
                alpha_ground: 0.8 as Float,
                frequency_ground: 4.0 as Float,
                alpha_excited: 1.2 as Float,
                frequency_excited: 5.0 as Float,
                ground_weight: vec![0.5 as TimeFloat, 0.5],
                excited_weight: vec![0.5 as TimeFloat, 0.5],
                distance: vec![0.25 as Float, 0.5],
            },
            filling: FillingTransitionInputs {
                characteristic_dose: vec![4.0 as Float, 4.0],
                dose_rate: 2.0 as Float,
                occupied_population: vec![0.2 as Float, 0.3],
                total_population: 1.0 as Float,
            },
            excited_energy_gap: vec![0.04 as Float, 0.05],
            s_frequency_e: vec![1.0e12, 1.0e9],
            s_frequency_g: vec![1.0e12, 1.0e9],
        };
        inputs
            .update_inputs(450.0, vec![0.4, 0.6])
            .expect("vector transition inputs should update");

        let rates = selected
            .calculate(&inputs)
            .expect("vector transition rates should calculate");

        assert!(rates.iter().all(|output| output.rate.len() == 2));
        assert!(rates
            .iter()
            .flat_map(|output| output.rate.iter())
            .all(|rate| rate.is_finite()));
    }

    #[test]
    fn runtime_selected_delocalised_equation_uses_one_common_interface() {
        let e_cb = vec![0.35 as Float, 0.45 as Float];
        let s_frequency = vec![1.0e12 as Float, 2.0e12 as Float];
        let temp = 450.0 as Float;

        let selected = DelocalisedRateEquationType::SecondOrder;
        let actual = selected
            .calculate(&e_cb, &s_frequency, &temp)
            .expect("runtime-selected equation should calculate rates");

        let expected = rate_equations::second_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
        )
        .expect("direct equation should calculate rates");

        assert_eq!(actual, expected);
    }

    #[test]
    fn delocalised_equation_can_be_selected_from_two_values() {
        assert_eq!(
            DelocalisedRateEquation::from_strs("ground", "first").unwrap(),
            DelocalisedRateEquation::Ground {
                re: DelocalisedRateEquationType::FirstOrder,
            },
        );

        assert_eq!(
            DelocalisedRateEquation::from_strs("excited", "second-order").unwrap(),
            DelocalisedRateEquation::Excited {
                re: DelocalisedRateEquationType::SecondOrder,
            },
        );

        assert_eq!(
            DelocalisedRateEquation::from_strs("both", "general:1.5").unwrap(),
            DelocalisedRateEquation::Both {
                re: DelocalisedRateEquationType::GeneralOrder { order: 1.5 },
            },
        );

        assert_eq!(
            DelocalisedRateEquation::from_strs("none", "first").unwrap(),
            DelocalisedRateEquation::None,
        );
    }

    #[test]
    fn delocalised_equation_can_be_selected_from_combined_text() {
        assert_eq!(
            "ground:first"
                .parse::<DelocalisedRateEquation>()
                .unwrap(),
            DelocalisedRateEquation::Ground {
                re: DelocalisedRateEquationType::FirstOrder,
            },
        );

        assert_eq!(
            "both:general:1.5"
                .parse::<DelocalisedRateEquation>()
                .unwrap(),
            DelocalisedRateEquation::Both {
                re: DelocalisedRateEquationType::GeneralOrder { order: 1.5 },
            },
        );

        assert_eq!(
            "none".parse::<DelocalisedRateEquation>().unwrap(),
            DelocalisedRateEquation::None,
        );
    }


    #[test]
    fn runtime_selection_accepts_vector_properties_and_scalar_temperature() {
        let e_cb: Vec<Float> = vec![0.35, 0.45];
        let s_frequency: Vec<Float> = vec![1.0e12, 2.0e12];
        let temp: Float = 450.0;

        let selected = DelocalisedRateEquationType::GeneralOrder { order: 1.5 };
        let actual = selected
            .calculate(&e_cb, &s_frequency, &temp)
            .expect("runtime selection should support scalar temperature broadcasting");

        let expected = rate_equations::general_order_delocalised_rate_equation(
            &e_cb,
            &s_frequency,
            &temp,
            &1.5,
        )
        .expect("direct equation should support scalar temperature broadcasting");

        assert_eq!(actual, expected);
    }


    #[test]
    fn nested_delocalised_transition_selects_ground_first_order() {
        let selected = DelocalisedRateEquation::Ground {
            re: DelocalisedRateEquationType::FirstOrder,
        };

        assert_eq!(
            selected,
            DelocalisedRateEquation::Ground {
                re: DelocalisedRateEquationType::FirstOrder,
            },
        );
    }

    #[test]
    fn nested_delocalised_transition_selects_both_general_order() {
        let selected = DelocalisedRateEquation::Both {
            re: DelocalisedRateEquationType::GeneralOrder { order: 1.5 },
        };

        assert_eq!(
            selected,
            DelocalisedRateEquation::Both {
                re: DelocalisedRateEquationType::GeneralOrder { order: 1.5 },
            },
        );
    }

    #[test]
    fn none_delocalised_equation_returns_scalar_zero_pair() {
        let selected = DelocalisedRateEquation::None;
        let inputs = DelocalisedTransitionInputs {
            e_cb_ground: 0.45_f32,
            frequency_ground: 1.0e12_f32,
            e_cb_excited: 0.55_f32,
            frequency_excited: 2.0e12_f32,
            temperature: 450.0,
            ground_weight: 1.0_f32,
            excited_weight: 1.0_f32,
        };

        let (ground, excited) = selected.calculate(&inputs);

        assert_eq!(ground.expect("ground zero should exist"), 0.0);
        assert_eq!(excited.expect("excited zero should exist"), 0.0);
    }

    #[test]
    fn none_delocalised_equation_returns_vector_zero_pair() {
        let selected = DelocalisedRateEquation::None;
        let e_cb_g: Vec<Float> = vec![0.35, 0.45, 0.55];
        let s_frequency_g: Vec<Float> = vec![1.0e12, 2.0e12, 3.0e12];
        let e_cb_e: Vec<Float> = vec![0.25, 0.35, 0.45];
        let s_frequency_e: Vec<Float> = vec![4.0e12, 5.0e12, 6.0e12];
        let temp: Float = 450.0;
        let ground_weight: Float = 1.0;
        let excited_weight: Float = 1.0;
        let inputs = DelocalisedTransitionInputs {
            e_cb_ground: e_cb_g,
            frequency_ground: s_frequency_g,
            e_cb_excited: e_cb_e,
            frequency_excited: s_frequency_e,
            temperature: temp,
            ground_weight,
            excited_weight,
        };

        let (ground, excited) = selected.calculate(&inputs);

        assert_eq!(
            ground.expect("ground zero vector should exist"),
            vec![0.0, 0.0, 0.0],
        );
        assert_eq!(
            excited.expect("excited zero vector should exist"),
            vec![0.0, 0.0, 0.0],
        );
    }


    #[test]
    fn filling_equation_can_be_selected_from_text() {
        assert_eq!(
            "basic".parse::<FillingRateEquation>().unwrap(),
            FillingRateEquation::Basic,
        );
        assert_eq!(
            "none".parse::<FillingRateEquation>().unwrap(),
            FillingRateEquation::None,
        );
        assert!("unknown".parse::<FillingRateEquation>().is_err());
    }

    #[test]
    fn basic_filling_selection_uses_filling_rate() {
        let d0: Float = 4.0;
        let d_dot: Float = 2.0;
        let n: Float = 0.25;
        let n_tot: Float = 1.0;
        let inputs = FillingTransitionInputs {
            characteristic_dose: d0,
            dose_rate: d_dot,
            occupied_population: n,
            total_population: n_tot,
        };

        let actual = FillingRateEquation::Basic
            .calculate(&inputs)
            .expect("basic filling equation should calculate a value");

        let expected = rate_equations::filling_rate(
            &d0,
            &d_dot,
            &n,
            &n_tot,
        )
        .expect("direct filling equation should calculate a value");

        assert_eq!(actual, expected);
    }

    #[test]
    fn none_filling_selection_returns_zero_matching_d0_shape() {
        let d0: Vec<Float> = vec![4.0, 5.0, 6.0];
        let d_dot: Float = 2.0;
        let n: Vec<Float> = vec![0.25, 0.5, 0.75];
        let n_tot: Float = 1.0;
        let inputs = FillingTransitionInputs {
            characteristic_dose: d0,
            dose_rate: d_dot,
            occupied_population: n,
            total_population: n_tot,
        };

        let actual = FillingRateEquation::None
            .calculate(&inputs)
            .expect("none filling equation should return zeros");

        assert_eq!(actual, vec![0.0, 0.0, 0.0]);
    }


    #[test]
    fn retrapping_selection_parses_named_components() {
        let selected = "cb_fi_gs_es"
            .parse::<RetrappingSelection>()
            .expect("retrapping selection should parse");

        assert_eq!(
            selected,
            RetrappingSelection {
                cb: true,
                filling: true,
                ground: true,
                excited: true,
            },
        );
        assert_eq!(selected.localised_rate_equation(), LocalisedRateEquation::Both);
    }

    #[test]
    fn transitions_from_strs_uses_named_retrapping_selection() {
        let selected = Transitions::from_strs(
            "cb_fi_gs",
            "both",
            "first",
            "ground",
            "basic",
        )
        .expect("valid transition configuration should parse");

        assert_eq!(
            selected,
            Transitions::FillCbRetrapping {
                transitions: TransitionsTypes{
                    delocalised: DelocalisedRateEquation::Both {
                    re: DelocalisedRateEquationType::FirstOrder,
                },
                   localised_recomb: LocalisedRateEquation::Ground,
                   localised_retrap: LocalisedRateEquation::Ground,
                   filling: FillingRateEquation::Basic,
                } 
            },
        );
    }

    #[test]
    fn cb_retrapping_requires_delocalised_equation() {
        let result = Transitions::from_strs(
            "cb",
            "none",
            "first",
            "none",
            "none",
        );

        assert!(result.is_err());
    }

    #[test]
    fn filling_retrapping_requires_filling_equation() {
        let result = Transitions::from_strs(
            "fi",
            "none",
            "first",
            "none",
            "none",
        );

        assert!(result.is_err());
    }

    #[test]
    fn retrapping_selection_accepts_compact_and_separated_forms() {
        let compact = "cbfigses"
            .parse::<RetrappingSelection>()
            .expect("compact retrapping selection should parse");
        let separated = "cb_fi_gs_es"
            .parse::<RetrappingSelection>()
            .expect("separated retrapping selection should parse");

        assert_eq!(compact, separated);
    }

    #[test]
    fn retrapping_selection_rejects_duplicate_and_unknown_components() {
        assert!("cb_cb".parse::<RetrappingSelection>().is_err());
        assert!("gsesgs".parse::<RetrappingSelection>().is_err());
        assert!("cb_unknown".parse::<RetrappingSelection>().is_err());
    }

    #[test]
    fn retrapping_selection_derives_all_localised_modes() {
        assert_eq!(
            "none".parse::<RetrappingSelection>().unwrap().localised_rate_equation(),
            LocalisedRateEquation::None,
        );
        assert_eq!(
            "gs".parse::<RetrappingSelection>().unwrap().localised_rate_equation(),
            LocalisedRateEquation::Ground,
        );
        assert_eq!(
            "es".parse::<RetrappingSelection>().unwrap().localised_rate_equation(),
            LocalisedRateEquation::Excited,
        );
        assert_eq!(
            "gses".parse::<RetrappingSelection>().unwrap().localised_rate_equation(),
            LocalisedRateEquation::Both,
        );
    }
}
