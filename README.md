# MCRustLum

MCRustLum is a Rust project for modelling charge-transport processes using
Monte Carlo simulations. It is primarily being developed for feldspar, but the
underlying approach can be applied to other crystal systems.

Analytical solvers are planned for a future release.

Detailed documentation can be found in the accompanying [website](https://dtuluminescencephysics.github.io/MCRustLUM/luminescence/)

## Project status

> [!IMPORTANT]
> MCRustLum is under active development

## Requirements

- A current [Rust toolchain](https://www.rust-lang.org/tools/install)
- Cargo, which is installed with Rust

## Building

Build every crate in the workspace with optimisations enabled:

```console
cargo build --workspace --release
```

## Testing

Run all unit and documentation tests:

```console
cargo test --workspace
```

## Input configuration

Simulation settings are grouped in [`input.toml`](input.toml). The application
copies this file into the experiment directory before starting a run.

## Workspace structure

| Crate | Purpose | Workspace Dependencies | External Dependencies |
| --- | --- | --- | --- |
| `common` | Crystal geometry, numerical operations, rate equations, and time-temperature profiles | — |`rand` 0.10.0, `rand_distr` 0.6.0, `ndarray` 0.17.2, `serde` 1.0 |
| `io` | Grouped simulation inputs, defaults, and TOML loading |`common` | `csv` 1.3, `flate2` 1.0, `plotters` 0.3.7, `serde` 1.0, `toml` 0.8 |
| `mc` | Monte Carlo simulation setup and state management | `common`, `io` | `rand` 0.10.0, `serde` 1.0, `rayon` 1.12.0 |
| `chronology` | Thermochronology simulation setup and state management | `common`, `io`, `mc` | `rand` 0.10.0, `statrs` 0.19, `serde` 1.0 |
| `luminescence` | Application entry point | `common`, `io`, `mc` `chronology` | — |

External dependency versions show the requirements declared in `Cargo.toml`; exact resolved versions are recorded in [`Cargo.lock`](Cargo.lock).
## Running

Run with an automatically numbered experiment directory:

```console
cargo run --release -p luminescence -- [Options]
```
### Options:
| Short Flag | Long Flag | Input | Purpose | Default |
| --- | ---| --- | --- | --- |
| -f | --folder |  NAME   | Experiment folder | automatic numbering |
| -r | --repetitions | NUMBER  | Repetitions | 1000 |
| -e | --experiments | NUMBER  | Experiments | 1 |
| -s | --seed | NUMBER  | Random seed | 0 |
| -t | --traps | NUMBER  | Minimum traps. | 800 |
| -mc | --montec | BOOLEAN | Run a Monte Carlo Simulation | TRUE |
| -tc | --thermochron | BOOLEAN | Run a Thermochronology experiment | FALSE |
| -p | --profiles | NUMBER  | Number of profiles to test during a thermochronology experiment | 10000 |
| -h | --help | — | Show help information | — |

The original single positional experiment-folder name is also supported."

For multiple threads: 
```console
RAYON_NUM_THREADS=#1 cargo run -p luminescence -- [Options]
```
## Output
The output behaviour needs to be overhauled to produce better outputs files 
that can be easily used in other graphing software or for further result 
interpretation. But the current behaviour is set as follows.
 
The output layout is `run/<name>/`, containing the copied `input.toml`, the
consolidated `average_fill.csv` and `average_event.csv`, six PNG plots, and a
`tmp/` directory for the per-repetition Monte Carlo files. The plots include
mean filling with standard-deviation limits, median filling with interquartile
limits, temperature against time, and aggregate event frequencies. Automatic
names use the first available `experiment_N`, starting with `experiment_1`.

## License

Copyright © 2026 Oliver Bramley, Technical University of Denmark.

MCRustLum is licensed under the
[GNU Affero General Public License version 3.0 only](LICENSE)
(`AGPL-3.0-only`).