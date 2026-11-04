// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

use crate::{profiles::{MoveKind, RandomWalkConstructor}};
use io::chrono_inputs::ChronologyInputs;
use common::numeric::Float;

pub struct AcceptanceWindow {
    pub move_kind: MoveKind,
    pub acceptance_range: (Float, Float),
    pub trials: usize,
    pub accepted: usize,
    pub stable_windows: usize,
}

impl AcceptanceWindow {
    pub fn new(move_kind: MoveKind, acceptance_range: (Float, Float)) -> Self {
        Self {
            move_kind,
            acceptance_range,
            trials: 0,
            accepted: 0,
            stable_windows: 0,
        }
    }

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

    pub fn stable_for(&self, required_windows: usize) -> bool {
        self.stable_windows >= required_windows
    }

    pub fn stable_windows(&self) -> usize {
        self.stable_windows
    }

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

pub struct BurnInApparatus  {
    pub acceptance_windows: Vec<AcceptanceWindow>,
    pub walker_alt: RandomWalkConstructor,
    pub window_rates: Vec<(Float, bool, usize,usize)>
}
impl BurnInApparatus {

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

    pub fn get_window(&mut self, move_kind: MoveKind) -> Result<&mut AcceptanceWindow, String>{
        self.acceptance_windows
            .iter_mut()
            .find(|trial_move| trial_move.move_kind == move_kind)
            .ok_or_else(|| format!("Missing acceptance window for {move_kind:?}"))

    }

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
