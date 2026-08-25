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

Simulation settings are grouped in [`input.toml`](input.toml). When that file
is not present, the program can construct the same input structures using its
built-in default values.

## Workspace structure

| Crate | Purpose |
| --- | --- |
| `common` | Crystal geometry, numerical operations, rate equations, and time-temperature profiles |
| `io` | Grouped simulation inputs, defaults, and TOML loading |
| `mc` | Monte Carlo simulation setup and state management |
| `luminescence` | Application entry point |

## Running

The application entry point is still being developed. A complete simulation
workflow and command-line interface will be documented here once available.

## License

Copyright © 2026 Oliver Bramley, Technical University of Denmark.

MCRustLum is licensed under the
[GNU Affero General Public License version 3.0 only](LICENSE)
(`AGPL-3.0-only`).
