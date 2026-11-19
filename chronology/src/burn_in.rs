// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Tune Bayesian walk proposal scales using per-move acceptance windows.
//!
//! Convergence requires consecutive observed windows inside every target range.
//! Birth and death corrections share scales; coordinate moves are tuned separately.

use crate::profiles::{MoveKind, RandomWalkConstructor};
use common::numeric::Float;
use io::chrono_inputs::ChronologyInputs;


/// Acceptance counts and consecutive stable-window count for one proposal kind.
pub struct AcceptanceWindow {
    /// Proposal kind whose acceptance is counted by this window.
    pub move_kind: MoveKind,
    /// Inclusive target acceptance-rate interval.
    pub acceptance_range: (Float, Float),
    /// Number of attempts observed in the current window.
    pub trials: usize,
    /// Number of successful attempts in the current window.
    pub accepted: usize,
    /// Number of consecutive observed windows within the target range.
    pub stable_windows: usize,
}

impl AcceptanceWindow {
    /// Start an empty observation window with an inclusive target acceptance range.
    pub fn new(move_kind: MoveKind, acceptance_range: (Float, Float)) -> Self {
        Self {
            move_kind,
            acceptance_range,
            trials: 0,
            accepted: 0,
            stable_windows: 0,
        }
    }

    /// Record one proposal attempt, incrementing the accepted count when successful.
    pub fn record(&mut self, accepted: bool) {
        self.trials += 1;
        self.accepted += usize::from(accepted);
    }

    /// Finish this move type's observations for one global window.
    /// A window with no attempts cannot count as stable.
    pub fn finish_window(
        &mut self,
    ) -> (Float, bool, usize, usize) {

        if self.trials == 0 {
            self.stable_windows = 0;
            return (0.0, false, 0, 0);
        }

        let trials = self.trials;
        let accepted = self.accepted;
        let rate = self.accepted as Float / self.trials as Float;
        let in_range = (self.acceptance_range.0..=self.acceptance_range.1).contains(&rate);
        self.stable_windows = if in_range { self.stable_windows + 1 } else { 0 };

        self.trials = 0;
        self.accepted = 0;

        (rate, in_range, trials, accepted)
    }

    /// Return whether at least `required_windows` consecutive windows met the target.
    pub fn stable_for(&self, required_windows: usize) -> bool {
        self.stable_windows >= required_windows
    }

    /// Return the current number of consecutive observed windows within the target range.
    pub fn stable_windows(&self) -> usize {
        self.stable_windows
    }

    /// Return target-midpoint minus observed acceptance when outside the target range.
    ///
    /// Returns zero inside the range and `None` when this move has no attempts.
    /// The correction is consumed before resetting the observation window.
    pub fn correction_factor(&self) -> Option<Float> {
        
        if self.trials == 0 {
            return None;
        }
        let rate = self.accepted as Float / self.trials as Float;
        let target = (self.acceptance_range.0 + self.acceptance_range.1) / 2.0;

        if rate < self.acceptance_range.0 || rate > self.acceptance_range.1 {
            Some(target - rate)
        } else {
            Some(0.0)
        }
    }



}

/// Adaptive proposal tuning state for Bayesian random-walk burn-in.
///
/// Keeps separate acceptance windows for birth, death, time, and temperature
/// moves, together with the walker whose scales are being tuned.
pub struct BurnInApparatus  {
    /// One observation window for each of the four proposal kinds.
    pub acceptance_windows: Vec<AcceptanceWindow>,
    /// Walker copy whose proposal scales are adapted during burn-in.
    pub walker_alt: RandomWalkConstructor,
    /// Finished-window tuples `(rate, in_range, trials, accepted)` per move.
    pub window_rates: Vec<(Float, bool, usize,usize)>
}
impl BurnInApparatus {

    /// Construct one acceptance window per move kind and an independently tunable walker.
    ///
    /// Uses the configured target ranges and proposal scales; walker construction
    /// errors are returned.
    pub fn new(inputs: &ChronologyInputs) -> Result<Self, String> {

        let acceptance_windows = vec![
            AcceptanceWindow::new(MoveKind::Birth, inputs.burn_in.birth_range),
            AcceptanceWindow::new(MoveKind::Death, inputs.burn_in.death_range),
            AcceptanceWindow::new(MoveKind::Time,  inputs.burn_in.time_range),
            AcceptanceWindow::new(MoveKind::Temp,  inputs.burn_in.temp_range)
        ];
        let walker_alt = RandomWalkConstructor::new(inputs)?;
        let window_rates = vec![(0.0, false, 0, 0);4];
        Ok(Self { 
            acceptance_windows,
            walker_alt,
            window_rates

        })

    }

    /// Borrow the observation window for a move kind, or report a missing window.
    pub fn get_window(&mut self, move_kind: MoveKind) -> Result<&mut AcceptanceWindow, String>{
        self.acceptance_windows
            .iter_mut()
            .find(|trial_move| trial_move.move_kind == move_kind)
            .ok_or_else(|| format!("Missing acceptance window for {move_kind:?}"))

    }

    /// Rescale proposal widths exponentially using the supplied acceptance corrections.
    ///
    /// Birth and death share their mean available correction. Time and temperature
    /// use their own corrections. Updated widths are clamped to `1e-8..=1e8`.
    pub fn adjust_walker(
        &mut self, 
        birth: Option<Float>, 
        death:Option<Float>, 
        time: Option<Float>, 
        temp: Option<Float>) 
    {
        let birth_death = match (birth, death) {
            (Some(birth), Some(death)) => Some((birth + death) / 2.0),
            (Some(value), None) | (None, Some(value)) => Some(value),
            (None, None) => None,
        };
        if let Some(birth_death) = birth_death {
            self.walker_alt.move_sigmas.time_birth = ( self.walker_alt.move_sigmas.time_birth 
                                                       * (-birth_death).exp()).clamp(1e-8, 1e8);
            self.walker_alt.move_sigmas.temp_birth = ( self.walker_alt.move_sigmas.temp_birth 
                                                   * (-birth_death).exp()).clamp(1e-8, 1e8);
        }
        if let Some(time) = time {
        self.walker_alt.move_sigmas.time = ( self.walker_alt.move_sigmas.time 
                                             * (-time).exp()).clamp(1e-8, 1e8);
        }

        if let Some(temp) = temp {
            self.walker_alt.move_sigmas.temp = ( self.walker_alt.move_sigmas.temp 
                                                 * (-temp).exp()).clamp(1e-8, 1e8);
        }

    }

/// Finish a global observation window and report whether tuning has converged.
///
/// Capture corrections before resetting counts, print per-move acceptance,
/// and require every move to meet `stable_window` consecutive target windows.
/// If not converged, adjust proposal scales. Missing move windows return errors.
/// `window_index` is zero-based and `maximum_no_windows` is used for progress output.
pub fn burn_in_operation(
    &mut self,
    window_index: usize,
    maximum_no_windows: usize,
    stable_window: usize,
) -> Result<bool, String> {

    // Capture corrections before finish_window resets the observation counts.
    let birth = self.get_window(MoveKind::Birth)?.correction_factor();
    let death = self.get_window(MoveKind::Death)?.correction_factor();
    let time = self.get_window(MoveKind::Time)?.correction_factor();
    let temp = self.get_window(MoveKind::Temp)?.correction_factor();

    for (w_rate, a_win) in self.window_rates.iter_mut().zip(self.acceptance_windows.iter_mut()){
        *w_rate = a_win.finish_window();
        if w_rate.2 > 0 {
              eprintln!(
                "Burn-in window {}/{}: {} acceptance {:.3} ({}/{}), target [{:.2}, {:.2}], stable windows {}",
                window_index + 1, maximum_no_windows, a_win.move_kind.to_name(),
                w_rate.0, w_rate.3, w_rate.2,
                a_win.acceptance_range.0, a_win.acceptance_range.1,
                if w_rate.1 {
                    a_win.stable_windows()
                } else {
                    0
                },
            );
        } else {
            eprintln!(
                "Burn-in window {}/{}: {} was not selected",
                window_index + 1, maximum_no_windows, a_win.move_kind.to_name(),
            );
        }
    }
   
    if self.acceptance_windows.iter().all(|a_win| a_win.stable_for(stable_window) ){
         eprintln!(
            "Burn-in converged after {} window(s). Tuned sigmas: time_birth={}, temp_birth={}, time={}, temp={}",
            window_index + 1,
            self.walker_alt.move_sigmas.time_birth, self.walker_alt.move_sigmas.temp_birth, 
            self.walker_alt.move_sigmas.time, self.walker_alt.move_sigmas.temp,
        );
        return Ok(true);
    }
    
    self.adjust_walker(birth, death, time, temp);
   
    Ok(false)
}

/// Print the tuned proposal scales and report reaching the configured window limit.
pub fn print_update(&self, window_index: usize, max_windows: usize) {
    if window_index + 1 == max_windows {
        eprintln!("Burn-in reached the maximum of {} window(s) without meeting the stability criterion.",max_windows);
    }
    eprintln!(
                "Current sigmas: time_birth={}, temp_birth={}, time={}, temp={}",
                self.walker_alt.move_sigmas.time_birth, self.walker_alt.move_sigmas.temp_birth,
                self.walker_alt.move_sigmas.time, self.walker_alt.move_sigmas.temp,
            );

}

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profiles::{MoveProbabilities, MoveSigmas};

    /// Construct a small walker with equal move weights and unit scales for burn-in tests.
    fn walker() -> RandomWalkConstructor {
        RandomWalkConstructor {
            age_range: (0.0, 1.0),
            temp_range: (0.0, 1.0),
            fixed_point_last: false,
            maintain_constraint_boxes: false,
            monotonicity: 'n',
            move_probs: MoveProbabilities::new(0.25, 0.25, 0.25, 0.25, 10),
            move_sigmas: MoveSigmas::new(1.0, 1.0, 1.0, 1.0),
        }
    }

    #[test]
    /// Verify that birth death corrections balance both moves.
    fn birth_death_corrections_balance_both_moves() {
        let kinds = [
            MoveKind::Birth,
            MoveKind::Death,
            MoveKind::Time,
            MoveKind::Temp,
        ];
        // Both below, both above, opposite sides, one in range, and no death trials.
        for (birth, death, expected_correction) in [
            (Some(5), Some(15), 0.15_f64),
            (Some(45), Some(35), -0.15),
            (Some(5), Some(35), 0.05),
            (Some(5), Some(45), 0.0),
            (Some(5), Some(25), 0.10),
            (Some(5), None, 0.20),
        ] {
            let mut burn = BurnInApparatus {
                acceptance_windows: kinds.iter().map(|&kind| AcceptanceWindow::new(kind, (0.2, 0.3))).collect(),
                walker_alt: walker(),
                window_rates: vec![(0.0, false, 0, 0); 4],
            };
            for (index, accepted) in [birth, death, Some(25), Some(25)].into_iter().enumerate() {
                if let Some(accepted) = accepted {
                    for trial in 0..100 {
                        burn.get_window(kinds[index]).unwrap().record(trial < accepted);
                    }
                }
            }
            assert!(!burn.burn_in_operation(0, 10, 2).unwrap());
            let walker = &burn.walker_alt;
            let expected = (-expected_correction).exp();
            assert!((walker.move_sigmas.time_birth - expected).abs() < 1e-12);
            assert!((walker.move_sigmas.temp_birth - expected).abs() < 1e-12);
            assert_eq!(walker.move_sigmas.time, 1.0);
            assert_eq!(walker.move_sigmas.temp, 1.0);
            // Finishing the window resets counts on the caller's windows.
            assert_eq!(burn.get_window(MoveKind::Birth).unwrap().finish_window().2, 0);
        }
    }

    #[test]
    /// Verify that finds moves in any order and reports convergence.
    fn finds_moves_in_any_order_and_reports_convergence() {
        let mut burn = BurnInApparatus {
            acceptance_windows: [MoveKind::Temp, MoveKind::Time, MoveKind::Death, MoveKind::Birth]
                .into_iter()
                .map(|kind| AcceptanceWindow::new(kind, (0.4, 0.6)))
                .collect(),
            walker_alt: walker(),
            window_rates: vec![(0.0, false, 0, 0); 4],
        };
        for index in 0..2 {
            for kind in [MoveKind::Birth, MoveKind::Death, MoveKind::Time, MoveKind::Temp] {
                let window = burn.get_window(kind).unwrap();
                assert_eq!(window.move_kind, kind);
                window.record(true);
                window.record(false);
            }
            assert_eq!(burn.burn_in_operation(index, 10, 2).unwrap(), index == 1);
        }
        assert_eq!(burn.walker_alt.move_sigmas.time_birth, 1.0);
        assert_eq!(burn.walker_alt.move_sigmas.temp_birth, 1.0);
        // An unobserved window breaks the consecutive stability requirement.
        assert!(!burn.burn_in_operation(2, 10, 2).unwrap());
        burn.acceptance_windows.retain(|window| window.move_kind != MoveKind::Birth);
        assert!(burn.get_window(MoveKind::Birth).is_err());
    }

    #[test]
    /// Verify that stability requires consecutive observed windows.
    fn stability_requires_consecutive_observed_windows() {
        let mut window = AcceptanceWindow::new(MoveKind::Birth, (0.4, 0.6));
        for _ in 0..2 {
            window.record(true);
            window.record(false);
            assert_eq!(window.finish_window(), (0.5, true, 2, 1));
        }
        assert!(window.stable_for(2));
        assert_eq!(window.finish_window(), (0.0, false, 0, 0));
        assert!(!window.stable_for(1));
    }
}
