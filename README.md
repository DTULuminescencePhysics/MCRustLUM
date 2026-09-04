# MCRustLum

MCRustLum is a Rust project for modelling charge-transport processes using
Monte Carlo simulations. It is primarily being developed for feldspar, but the
underlying approach can be applied to other crystal systems.

Analytical solvers are planned for a future release.

Detailed documentation can be found in the accompanying [website](https://dtuluminescencephysics.github.io/MCRustLUM/luminescence/)

## Project status

> [!IMPORTANT]
> MCRustLum is under active development and does not yet provide a complete
> runnable simulation.

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
cargo run -p luminescence
```

Or provide a folder name:

```console
cargo run -p luminescence -- my_experiment
```

The output layout is `run/<name>/`, containing the copied `input.toml`, the
consolidated `average_fill.csv`, and a `tmp/` directory for the per-repetition
Monte Carlo files. Automatic names use the first available `experiment_N`,
starting with `experiment_1`.

## License

Copyright © 2026 Oliver Bramley, Technical University of Denmark.

MCRustLum is licensed under the
[GNU Affero General Public License version 3.0 only](LICENSE)
(`AGPL-3.0-only`).
