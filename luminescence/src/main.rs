// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Command-line entry point for a luminescence Monte Carlo run.
//!
//! The program creates a directory under `run/`, copies `input.toml` into it,
//! and runs the simulation from that directory.

use std::error::Error;
use std::ffi::OsString;
use std::time::Instant;

const HELP: &str = "Usage: luminescence [OPTIONS]

Options:
  -f, --folder      NAME    Experiment folder (default: automatic numbering)
  -r, --repetitions NUMBER  Repetitions (default: 1000)
  -e, --experiments NUMBER  Experiments (default: 1)
  -s, --seed        NUMBER  Random seed (default: 0)
  -t, --traps       NUMBER  Minimum traps (default: 800)
  -mc               BOOLEAN Run a Monte Carlo Simulation (default: true)
  -tc               BOOLEAN Run a Thermochronology experiment
  -p --profiles     NUMBER  Number of profiles to test during a thermochronology experiment (default: 10000) 
  -h, --help                Show this help

The original single positional experiment-folder name is also supported.";

/// Prepare the requested run directory and execute the default simulation workflow.
fn main() -> Result<(), Box<dyn Error>> {    
    let Some(arguments) = parse_arguments(std::env::args_os().skip(1))? else {
        println!("{HELP}");
        return Ok(());
    };
    io::filesystem::prepare_experiment_directory(arguments.folder.as_deref())?;
    common::random::set_seed(arguments.seed);

    let total_start = Instant::now();

    if arguments.mc{
        let start = Instant::now();
        mc::monte_carlo_run(
            arguments.repetitions, arguments.experiments, arguments.traps, None,
        )?;
        eprintln!("Monte Carlo Complete:             {:?}", start.elapsed());
    }
    if arguments.tc{
        let start = Instant::now();
        chronology::thermo_chronology_experiment(
            arguments.repetitions, arguments.experiments, arguments.traps, arguments.profiles,true)?;
        eprintln!("Thermochronology Complete:             {:?}", start.elapsed());
    }
    
    eprintln!("total:             {:?}", total_start.elapsed());

    Ok(())

}


#[derive(Debug)]
struct RunArguments {
    folder: Option<OsString>,
    repetitions: usize,
    experiments: usize,
    seed: u64,
    traps: usize,
    mc: bool,
    tc: bool,
    profiles:usize
}

/// Return `None` for help, before creating any experiment directory.
fn parse_arguments(arguments: impl IntoIterator<Item = OsString>) -> Result<Option<RunArguments>, String> {
    let mut arguments = arguments.into_iter();
    let mut parsed = RunArguments {
        folder: None,
        repetitions: 1000,
        experiments: 1,
        seed: 0,
        traps: 800,
        mc: true,
        tc: false,
        profiles:10000,
    };

    while let Some(argument) = arguments.next() {
        let flag = argument.to_str().unwrap_or("");
        match flag {
            "-h" | "--help" => return Ok(None),
            "-f" | "--folder" | "-r" | "--repetitions" | "-e" | "--experiments" |
            "-s" | "--seed"   | "-t" | "--traps" | "-mc" |"--montec"  |
            "--thermochron"|"-tc" | "-p" | "--profiles" => {
                let value = arguments.next().ok_or_else(|| format!("missing value for {flag}"))?;
                if value.to_string_lossy().starts_with('-') {
                    return Err(format!("expected a value for {flag}, got {value:?}"));
                }
                match flag {
                    "-f" | "--folder" => {
                        if parsed.folder.is_some() {
                            return Err("experiment folder was specified more than once".into());
                        }
                        parsed.folder = Some(value);
                    }
                    "-r" | "--repetitions" => parsed.repetitions = parse_number(flag, value)?,
                    "-e" | "--experiments" => parsed.experiments = parse_number(flag, value)?,
                    "-s" | "--seed" => parsed.seed = parse_number(flag, value)?,
                    "-t" | "--traps" => parsed.traps = parse_number(flag, value)?,
                    "-p" | "--profiles" => parsed.profiles = parse_number(flag, value)?,
                    "-mc"| "--montec" => parsed.mc = parse_number(flag, value)?,
                    "-tc"| "--thermochron" => parsed.tc = parse_number(flag, value)?,
                   
                    _ => unreachable!(),
                }
            }
            _ if !argument.to_string_lossy().starts_with('-') && parsed.folder.is_none() => {
                parsed.folder = Some(argument);
            }
            _ => return Err(format!("unknown argument {argument:?}; use --help for usage")),
        }
    }

    if parsed.repetitions == 0 || parsed.experiments == 0 || parsed.traps == 0 {
        return Err("repetitions, experiments, and traps must be positive".into());
    }
    Ok(Some(parsed))
}

fn parse_number<T: std::str::FromStr>(flag: &str, value: OsString) -> Result<T, String> {
    value.to_str().and_then(|text| text.parse().ok())
        .ok_or_else(|| format!("invalid integer for {flag}: {value:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(arguments: &[&str]) -> Result<Option<RunArguments>, String> {
        parse_arguments(arguments.iter().map(OsString::from))
    }

    #[test]
    fn parses_flags_in_any_order() {
        let args = parse(&["-s", "42", "-t", "900", "-e", "3", "-f", "my run", "-r", "2000"])
            .unwrap().unwrap();
        assert_eq!(args.folder, Some(OsString::from("my run")));
        assert_eq!((args.repetitions, args.experiments, args.seed, args.traps), (2000, 3, 42, 900));
        let args = parse(&["--folder", "long", "--repetitions", "2", "--experiments", "4",
            "--seed", "5", "--traps", "6"]).unwrap().unwrap();
        assert_eq!((args.repetitions, args.experiments, args.seed, args.traps), (2, 4, 5, 6));
    }

    #[test]
    fn preserves_defaults_and_supports_help_and_original_folder_argument() {
        let args = parse(&[]).unwrap().unwrap();
        assert!(args.folder.is_none());
        assert_eq!((args.repetitions, args.experiments, args.seed, args.traps), (1000, 1, 0, 800));
        assert_eq!(parse(&["old_folder"]).unwrap().unwrap().folder, Some("old_folder".into()));
        assert!(parse(&["--help"]).unwrap().is_none());
        assert!(parse(&["-h"]).unwrap().is_none());
    }

    #[test]
    fn rejects_missing_invalid_and_unknown_arguments() {
        for args in [vec!["-r"], vec!["-r", "abc"], vec!["-e", "0"], vec!["-t", "-1"],
            vec!["-s", "-1"], vec!["-r", "0"], vec!["-t", "0"], vec!["--unknown"],
            vec!["-f", "-r", "2"], vec!["one", "two"], vec!["one", "-f", "two"]] {
            assert!(parse(&args).is_err(), "should reject {args:?}");
        }
    }
}
