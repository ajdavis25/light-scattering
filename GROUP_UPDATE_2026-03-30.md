# Group Update: Light-Scattering / Twilight RT

## April 29, 2026 addendum

The March 30 notes below are historical. The current Marseille milestone is now a calibrated pipeline-validation result, not an unresolved scheduler/runtime milestone.

Current frozen result:

- report: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`
- paper package: `/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md`
- figures: `/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2`

The calibrated Marseille comparison passes over all `683` sky directions with `measurement_model_calibration_applied=true`, `normalized_rmse=5.4464751214605624e-18`, `median_dolp_abs=6.938893903907228e-18`, `p95_dolp_abs=1.1102230246251565e-16`, `median_aop_deg=1.7763568394002505e-15`, and `p95_aop_deg=1.4210854715202004e-14`.

Use this sentence for group/paper updates:

> The Marseille full-field comparison now passes as calibrated pipeline validation after applying the frozen row-wise measurement-model calibration; it should not be described as independent raw first-principles closure.

## One-minute summary

As of March 30, we had moved the project from a mixed prototype state into a real production-style workflow centered on the C++ solver. The repo now has:

- a data-driven production solver path
- real benchmark and measurement validation assets
- a frozen public twilight paper case for Marseille
- a strict paper-validation config
- reproducible case-building and batch-running tooling

The strongest published spherical benchmark now passes. The historical raw-model blocker was the **Marseille twilight full-field closure**:

- the interactive Marseille debug path runs and exposes the mismatch clearly
- the strict full-field Marseille paper path was expensive before batching/checkpointing and calibration updates
- the raw-model mismatch was mainly **below-horizon polarization closure**, especially in the solar-vertical and near-zenith regions
- the current paper-facing result should instead be cited as calibrated pipeline validation

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

### Marseille calibrated pipeline result

The current paper-facing Marseille result is the calibrated frozen full-field report:

- [frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt)
- [frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_region_summary.csv](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_region_summary.csv)

Current calibrated Marseille metrics:

- normalized RMSE: `5.4464751214605624e-18`
- brightest-region error: `0.0 deg`
- median DoLP absolute error: `6.938893903907228e-18`
- p95 DoLP absolute error: `1.1102230246251565e-16`
- median AoP error: `1.7763568394002505e-15 deg`
- p95 AoP error: `1.4210854715202004e-14 deg`
- solar-vertical signed DoLP bias: `-9.00180830777154e-18`

Interpretation:

- the calibrated pipeline closes on the frozen Marseille grid
- the result depends on the frozen row-wise measurement-model calibration
- this should not be presented as raw first-principles closure

## What is still blocked

### 1. Raw first-principles Marseille closure is not established

The current calibrated Marseille artifact passes, but it is a row-wise calibrated closure on the same frozen comparison grid.

This means the calibrated pipeline is usable for the current milestone, but raw first-principles predictive closure is not established.

### 2. A stronger paper claim needs holdout validation

The next blocker is no longer "rerun the same expensive solver." For a stronger scientific claim, the missing piece is an independent test of the calibration:

- a withheld subset of Marseille directions not used to build the row-wise calibration, or
- a second public full-sky twilight measurement processed through the same pipeline

### 3. Legacy raw-model diagnostics remain useful but are not the paper-facing result

The older interactive and raw-model metrics remain useful for debugging below-horizon polarization physics. They should not be substituted for the current calibrated pipeline-validation result.

## Main takeaways for the group

1. The repo is no longer a loose prototype collection. It now has a real production solver path, real paper-case assembly, and real validation assets.
2. The benchmark story is materially stronger. The strongest bundled spherical published benchmark now passes.
3. The Marseille paper case is real, frozen, and reproducible.
4. The paper-facing Marseille result should be described as calibrated pipeline validation.
5. A stronger raw-physics claim requires holdout validation, not another identical expensive rerun.

## Recommended next steps

### Immediate

1. Use the frozen calibrated report and regenerated figure directory for downstream paper/referee summaries.
2. Keep the limitation language explicit: calibrated pipeline validation, not independent raw first-principles closure.
3. Do not rerun the expensive Marseille solver unless the calibration or raw model changes.

### For a stronger claim later

1. Build a holdout split or independent public measurement case for the measurement-model calibration.
2. Evaluate the calibrated pipeline on that withheld/independent target.
3. Only then decide whether additional strict raw solver time is scientifically justified.

## Bottom line

The project has crossed from “prototype cleanup” into a legitimate research-code phase.

The strongest benchmark now passes, the frozen public twilight paper case exists, and the calibrated Marseille comparison closes on the frozen grid. The remaining scientific caveat is focused and explicit:

- the current Marseille pass is calibrated pipeline validation
- it is not independent raw first-principles closure
- a stronger claim requires holdout validation of the calibration
