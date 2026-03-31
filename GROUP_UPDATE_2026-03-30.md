# Group Update: Light-Scattering / Twilight RT

## One-minute summary

We have moved the project from a mixed prototype state into a real production-style workflow centered on the C++ solver. The repo now has:

- a data-driven production solver path
- real benchmark and measurement validation assets
- a frozen public twilight paper case for Marseille
- a strict paper-validation config
- reproducible case-building and batch-running tooling

The strongest published spherical benchmark now passes. The main remaining blocker is not infrastructure anymore. It is the **Marseille twilight full-field closure**:

- the interactive Marseille debug path runs and exposes the mismatch clearly
- the strict full-field Marseille paper path is still too expensive to complete as currently configured
- the current mismatch is still mainly **below-horizon polarization closure**, especially in the solar-vertical and near-zenith regions

## What changed

### 1. The C++ solver is now the scientific source of truth

The active solver path is now:

- [monte_carlo_cpp/src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp)
- [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp)

Key changes:

- production interface is now `SimulationConfig -> SkyResult`
- observer-centered fisheye bins are the native output
- order-resolved outputs now exist:
  - first order
  - second order
  - third-plus remainder
  - second-order `RR / AR / RA / AA` breakdown

### 2. Inputs are now data-driven instead of placeholder-only

The solver now reads:

- atmosphere profiles from CSV
- solar spectrum tables
- aerosol phase / Mueller lookup tables
- gas absorption tables for `O3`, `O2`, `O4`, `H2O`, and `NO2`
- instrument response curves

Key files:

- [monte_carlo_cpp/src/Atmosphere.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.cpp)
- [monte_carlo_cpp/src/WavelengthHandling.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.cpp)
- [monte_carlo_cpp/src/PhaseFunctions.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.cpp)
- [monte_carlo_cpp/src/Polarization.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Polarization.cpp)
- [monte_carlo_cpp/src/SurfaceReflection.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/SurfaceReflection.cpp)

### 3. Paper-path case assembly now exists

We now have a real frozen public twilight case:

- [monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg)
- [monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg)
- [monte_carlo_cpp/config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg)

The frozen case package is under:

- [monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z)

Builder:

- [monte_carlo_cpp/tools/build_paper_case_inputs.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/build_paper_case_inputs.py)

Current Marseille case assembly includes:

- Open-Meteo pressure-level thermodynamics
- Open-Meteo air-quality fields
- AERONET direct-sun AOD
- AERONET inversion-backed aerosol fallback
- public-doc-constrained IMX250MYR blue-channel response proxy

From the current provenance:

- timestamp: `2022-08-15T19:14:13Z`
- solar zenith: `96.0071 deg`
- solar azimuth: `295.9851 deg`
- AERONET inversion usage: `same_day_fallback`
- selected inversion time: `2022-08-15T15:53:27Z`
- inversion quality level: `lev15`

Relevant file:

- [paper_case_provenance.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/paper_case_provenance.json)

### 4. Validation is materially stronger now

We now have:

- executable unit/physics checks
- scalar and vector sanity benchmarks
- published spherical-vector benchmark support
- measurement-case evaluation with pointwise and regionwise diagnostics

The strongest spherical published benchmark currently passes:

- [benchmark_zawada_spherical_vector_multiple_smoke.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/benchmark_zawada_spherical_vector_multiple_smoke.txt)

Current passing Zawada multiple-scattering benchmark metrics:

- median intensity relative error: `0.00457567`
- p95 intensity relative error: `0.012131`
- median DoLP absolute error: `0.0199496`
- p95 DoLP absolute error: `0.0353911`

This is an important milestone: the solver is no longer failing the strongest bundled spherical benchmark.

### 5. Marseille measurement reduction and diagnostics now exist

The Marseille measurement is no longer just a candidate dataset. It has been reduced into the validator format:

- [measurement_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference.csv)
- [measurement_reduction.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reduction.json)

We also now have:

- strict subset reference for fast physics iteration
- AoP-aware comparison output
- regionwise Marseille summaries
- order-fraction diagnostics

Key files:

- [monte_carlo_cpp/tools/reduce_marseille_case.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/reduce_marseille_case.py)
- [monte_carlo_cpp/src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp)

## What works now

### Benchmarks

- Published spherical Zawada benchmark passes
- Unit tests and executable physics tests exist
- Paper-path input assembly is real and reproducible

### Marseille debug path

The interactive Marseille full-field debug path runs and produces useful diagnostics:

- [frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt)
- [frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv)

Current interactive Marseille metrics:

- normalized RMSE: `0.155497`
- brightest-region error: `23.2398 deg`
- median DoLP absolute error: `0.169169`
- p95 DoLP absolute error: `0.567559`
- median AoP error: `40.6645 deg`
- p95 AoP error: `83.4536 deg`
- solar-vertical signed DoLP bias: `0.374188`

Interpretation:

- the solver no longer collapses to zero field
- the earlier discrete AoP quadrant/sign issue was materially reduced
- second-order and higher-order contributions are now present and measurable
- the remaining mismatch is still too large for paper use

## What is still blocked

### 1. The strict Marseille paper gate is not tractable yet

The full strict Marseille measurement batch did **not** finish. After running for multiple days, it still had not exited the Marseille measurement phase, and `paper_validation.cfg` never started.

This was not a launcher bug. It was a runtime-budget problem in the strict measurement config.

### 2. The strict paper config is orders of magnitude heavier than the interactive debug surface

For the same `683` Marseille directions:

- spectral bands:
  - strict: `46`
  - interactive: `3`
- solar disk nodes:
  - strict: `7`
  - interactive: `1`
- photons per bin:
  - strict: `256`
  - interactive: `1`

Approximate work increase relative to the interactive full-field debug config:

- Monte Carlo remainder: about `3925x`
- deterministic second-scatter outer quadrature: about `34x`
- deterministic second-scatter including ray depth and solar-disk sampling: about `358x`

That is the main runtime diagnosis.

### 3. Paper-level Marseille closure is still not achieved

The current main science blocker is still:

- below-horizon polarization closure
- especially the solar-vertical mid-zenith and near-zenith regions

The solver still over-polarizes and under-fills key parts of the Marseille field, even though the infrastructure and diagnostics are now much better.

## Main takeaways for the group

1. The repo is no longer a loose prototype collection. It now has a real production solver path, real paper-case assembly, and real validation assets.
2. The benchmark story is materially stronger. The strongest bundled spherical published benchmark now passes.
3. The Marseille paper case is real, frozen, and reproducible.
4. The remaining blocker is no longer “missing infrastructure.” It is:
   - strict Marseille runtime tractability
   - Marseille twilight physics closure
5. We are now at the stage where solver/runtime optimization and physics closure matter more than repo scaffolding.

## Recommended next steps

### Immediate

1. Add progress and timing instrumentation to the Marseille measurement runner.
2. Add a strict-paper profiling subset that preserves paper physics but reduces the number of scored directions.
3. Optimize the strict Marseille path before relaunching the full paper gate.

### After runtime is under control

1. Re-run the strict Marseille full-field measurement case.
2. Then run [paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg).
3. Use the Marseille region diagnostics to target the remaining solar-vertical and near-zenith polarization errors.

## Bottom line

The project has crossed from “prototype cleanup” into a legitimate research-code phase.

The strongest benchmark now passes, the frozen public twilight paper case exists, and the measurement mismatch is now well exposed. The remaining work is focused and technical:

- make the strict Marseille run tractable
- close the Marseille twilight polarization mismatch
- then promote the paper gate from infrastructure-complete to scientifically passing
