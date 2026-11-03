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

| Crate | Purpose |
| --- | --- |
| `common` | Crystal geometry, numerical operations, rate equations, and time-temperature profiles |
| `io` | Grouped simulation inputs, defaults, and TOML loading |
| `mc` | Monte Carlo simulation setup and state management |
| `luminescence` | Application entry point |

## Running

Run with an automatically numbered experiment directory:

```console
cargo run --release -p luminescence -- [Options]
```
Options:
  -f, --folder      NAME    Experiment folder (default: automatic numbering)
  -r, --repetitions NUMBER  Repetitions (default: 1000)
  -e, --experiments NUMBER  Experiments (default: 1)
  -s, --seed        NUMBER  Random seed (default: 0)
  -t, --traps       NUMBER  Minimum traps (default: 800)
  -mc               BOOLEAN Run a Monte Carlo Simulation (default: true)
  -tc               BOOLEAN Run a Thermochronology experiment  
  -p, --profiles    NUMBER  Number of profiles to test during a thermochronology experiment (default: 10000)
  -h, --help                Show this help

The original single positional experiment-folder name is also supported."

For multiple threads: 
```console
RAYON_NUM_THREADS=#1 cargo run -p luminescence -- [Options]
```


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