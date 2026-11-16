# Understanding MCRustLum

MCRustLum follows trapped electrons through a crystal while a prescribed
thermal history changes transition rates. The Monte Carlo workflow predicts
filling and event histories. The chronology workflow proposes thermal histories
and compares their predicted final fills and ages with measurements.

Detailed documentation can be found in the accompanying [website](https://dtuluminescencephysics.github.io/MCRustLUM/luminescence/).

## Where to start reading

| Question | Entry point |
| --- | --- |
| What happens when I run the application? | [`main`](../luminescence/src/main.rs) parses arguments, prepares the experiment directory, seeds randomness, and selects workflows. |
| How are configuration names and units resolved? | [`io::inputs`](../io/src/inputs.rs) loads named entities and produces `ResolvedExperiment` values. |
| How is a Monte Carlo run assembled? | [`MonteCarloSimulation`](../mc/src/system_setup.rs) builds geometry, profiles, trap generators, occupancy, and transitions. |
| How does one repetition execute? | [`MCExperiment`](../mc/src/experiment.rs) constructs independent mutable state and starts the kinetic runner. |
| How is the next physical event chosen? | [`calculate_times`](../mc/src/calculate_times.rs) builds candidates, samples lifetimes, and applies the earliest event. |
| Where are the physical equations? | [`rate_equations`](../common/src/rate_equations.rs) implements the equations; [`rate_equation_selection`](../common/src/rate_equation_selection.rs) selects enabled models. |
| How are thermal histories inferred? | [`ChronologyExperiment`](../chronology/src/system_setup.rs) connects proposals, forward simulation, acceptance, and optional burn-in. |
| How are outputs interpreted? | [`mc::average`](../mc/src/average.rs) consolidates trajectories; [`chronology::analyse`](../chronology/src/analyse.rs) builds accepted-profile occupancy grids. |

## Monte Carlo workflow

1. `io::inputs::read_inputs` reads the requested experiments in
   `shared_settings.experiment_list` order. It resolves named trap, cube, band,
   filling, transition, and time-temperature definitions. Experiment-specific
   references override the shared defaults.
2. The loader normalizes rate parameters to each profile's time unit, spatial
   parameters to the cube height's metric prefix, and profile temperatures to
   Kelvin. Quantities can retain optional distributions for later sampling.
3. `MonteCarloSimulation::new` builds reusable configuration. Geometry uses trap
   densities to derive site counts and can enlarge the crystal to reach the
   requested minimum trap population. Geometry inputs `hole_count` and
   `bandtail_count` are passed as counts **per trap** to this density-based path.
4. Each `MCExperiment` receives independently randomized site positions,
   occupancy partitions, profile state, and a reproducible random stream.
   `TrapParameterLayout` maps site identifiers to uniform, individual, or shared
   family parameter records.
5. The kinetic runner rebuilds enabled event candidates at the current state.
   Competing exponential waiting times select the next event. A profile-step
   candidate limits travel to the next temperature/control-point boundary so
   rates can be recalculated as the temperature changes.
6. Applying an event updates trap and hole occupancy. Conduction-band release
   requires a further competition between recombination and retrapping
   destinations. An aggregate filling event resolves its selected sites when
   applied. The standard runner writes state/event records in compressed
   batches; the final-ratio runner returns only the endpoint fill.
7. `mc::average` merges repetition times, carries filling states forward, and
   produces filling statistics, cumulative repetition statistics, and event
   bins. `io::plotting` renders the consolidated results.

`PlaceAvailability` provides the occupancy bookkeeping: available trap IDs are
electron occupied, while available hole IDs contain holes that can recombine.
It stores a permutation and reverse lookup so changing a known site's
availability does not require scanning the population.

## Chronology workflow

1. `io::chrono_inputs::read_chrono_inputs` reads `chron_input.toml`, filling
   omitted fields from defaults. Construction and execution perform further
   validation; loading alone does not validate a complete inference setup.
2. `RandomProfileConstructor` samples constraint-box endpoints and recursively
   inserts intermediate nodes. `ProfileType::Walk` instead modifies the current
   profile using birth, death, time, or temperature proposals after initialization.
3. Proposal validation checks structure, finite coordinates, bounds, strict
   spacing, and temperature monotonicity. A rejected constraint proposal consumes
   a trial without running the forward model. Profile generation codes are: 
      `R` - Random value generation each profile is independent 
      `W` - Walk through time/temperature space each profile is based on the previous.
   Monotonicity codes are: 
      `I` - strictly increasing temperature values, 
      `D` - strictly Decreasing temperature values, 
      `V` - temperature values can vary.
   If a random profile generation is selected the number of halving segments inserted
   between each constraint box is set by `intermediate_node_halving_num`. The number of
   intermediate points is calculated n^2-1 so 3 halving points would give 8 intermediate 
   points. Each point is then calculated with one of three regimes:
      `E` - Episodic, 1 random number, 
      `I` - Intermediate, average of 2 numbers, 
      `G` - Gradual average of 4 numbers.
4. The forward Monte Carlo model returns final fill ratios indexed by experiment
   and then repetition. `AcceptanceConstructor` compares their means with the
   corresponding measured fills, using squared uncertainty-normalized residuals.
   The age comparison currently uses the first measured age and its uncertainty.
5. Frequentist acceptance combines individual p-values with Fisher's method and
   rejects below 0.05. Independent Bayesian acceptance uses the trial log
   likelihood. Bayesian walk acceptance compares current and trial likelihoods
   and includes move-specific corrections; accepted moves update the current
   profile. The stored Bayesian `log_posterior` field remains a zero placeholder.
6. Optional Bayesian walk burn-in monitors each move separately. Consecutive
   target-range windows establish convergence; exponential corrections adapt
   proposal scales. Production then buffers accepted profiles in batches of
   1,000 and flushes the remainder at completion.
7. Chronology output stores units followed by accepted-profile batches in
   `tmp/chronology_results.bin.gz`. `analyse` clips profile segments to grid
   bounds and counts each profile at most once per crossed cell. The resulting
   CSV and plots show occupancy, representative paths, and quantile bands.

## Units and conventions

- `TimeTemperature` retains numerical time values in their supplied profile
  unit. `Ka` and `Ma` select decreasing geological age; other units select
  increasing elapsed time. Rates must use the matching time unit. The explicit
  `convert_to_seconds` helper is available when seconds are required elsewhere.
- The low-level profile type preserves temperature numbers; it does not itself
  convert Celsius to Kelvin. The simulation input loader converts its profile
  temperatures. Directly constructed forward-model profiles must supply Kelvin.
- Chronology records carry paired vectors without embedded unit metadata. Their
  configured/stored units must agree with the supplied values, measurements,
  and overlays. Changing a setup unit label does not convert those vectors.
- `TimeUnit::convert_unit` treats `self` as the target. Temperature conversion
  dispatches on the `(self, unit)` pair. Read their method documentation before
  using these low-level conversion interfaces.
- The metric prefix parser falls back to `STANDARD` for unknown names. The
  current `GIGA` volume multiplier returns 1.0 to stop code failing but 
  should not be used. These behaviors are documented where implemented; 
  this documentation update does not change calculations.

## Read and maintain the generated reference

Build documentation for all crates, including private implementation details:

```console
cargo doc --workspace --no-deps --document-private-items
```

Open `target/doc/luminescence/index.html` for the application entry point,
`target/doc/common/index.html` for the physical model, or
`target/doc/chronology/index.html` for inference. Every documented item has a
source link so readers can move between explanation and implementation.

Run the source coverage check and executable documentation examples:

```console
python3 docs/check_documentation.py
cargo test --workspace --doc
```

The source check includes private/test declarations and recognizes the
workspace's macro documentation forwarding. It is a lightweight source check,
not a general Rust parser. Rustdoc independently validates expanded public
items, fields, variants, and links. The documentation workflow runs both checks,
builds private-item documentation, and tests the examples before publishing.
