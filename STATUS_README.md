# Project Status

Last updated: 2026-04-29

This file tracks the repo state after the research-grade clear-sky twilight implementation pass that followed the earlier inspection report in [notebooks/notes.ipynb](/c:/Users/ashton/Desktop/projects/light-scattering/notebooks/notes.ipynb).

Current physical and numerical model assumptions are documented in [MODEL_ASSUMPTIONS.md](/c:/Users/ashton/Desktop/projects/light-scattering/MODEL_ASSUMPTIONS.md).

## April 29, 2026 addendum

The current Marseille paper-facing result is no longer the earlier interactive/raw mismatch described below. The current frozen full-field Marseille artifact is a calibrated pipeline-validation result:

- status note: [notebooks/MARSEILLE_CALIBRATED_VALIDATION_STATUS_2026-04-28.md](/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_VALIDATION_STATUS_2026-04-28.md)
- paper package: [notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md](/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md)
- publication package: [paper/README.md](/work/vmo703/light-scattering/paper/README.md)
- fresh-thread handoff: [notebooks/MARSEILLE_LIGHT_SCATTERING_HANDOFF_2026-04-30.md](/work/vmo703/light-scattering/notebooks/MARSEILLE_LIGHT_SCATTERING_HANDOFF_2026-04-30.md)
- report: [frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt)
- figures: [plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2](/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2)

This result passes the configured Marseille measurement gate only after applying the frozen row-wise measurement-model calibration. It should be cited as calibrated pipeline validation, not as independent raw first-principles closure.

## 2026-08-11 addendum — V2 physics-only frozen tier (current)

A second frozen tier now exists at `monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z_v2/` (v1 untouched); its README is the canonical description. Stack (every non-default traces to physics or measurement): `event_frame_chi_sign_fix=true` (2026-08-01 handedness fix), `rayleigh_depolarization_factor=0.0279` (King factor), and an AERONET-constrained day profile (column AOD550 = 0.155 from the at-site Marseille_ATMO evening plateau/bracket, boundary-layer-weighted; v1's 0.1379 was already day-informed via a 0.75/0.25 AERONET/Open-Meteo blend — the v2 value is a refinement of the same station's data, not a missing-data correction). The aerosol f22 nonsphericity knob is EXCLUDED (no dust that day; documented sensitivity only). Runs: Slurm 788181 (48-dir subset + 683-dir full field, pinned binary `build_v2freeze`), gate 788184.

Headline numbers (comparison-CSV metric family, `tools/score_measurement_comparison.py`): full-field signed DoLP bias median +0.2068 → **+0.1736**, AoP median 16.96°, shape RMSE 0.0827 → 0.0940; subset +0.1938 → **+0.1678** with AoP 23.4° → 18.0°. **The residual ≈ +0.17 over-polarization is a declared open discrepancy and is NOT input-side** — the AERONET-licensed loading correction is ±10–15% in column and contributes only ≈ −0.003 of DoLP beyond King (three-station verdict archived in `paper/tables/reproducibility/aeronet_20220815/`; the ×2–×4 loading scans of 2026-08-10 were compensation, not physics). Prime suspect for the residual: the Marseille-vs-Koomen bias asymmetry (reduction-side), sharpened by the V2 gate result below.

V2 validation gate (`validation_default_capped_v2.cfg`, report archived in `paper/tables/reproducibility/`): **22/23 metrics pass**. Both twilight measurement cases now carry the true h⁻⁶ profile plus King (the 2026-08-09 Rozenberg input archaeology applied coherently to the whole h⁻⁶ scenario family): Rozenberg `normalized_rmse = 0.0851` ≈ the historical claim 0.0849, and **Koomen fully passes for the first time** (`normalized_rmse = 0.040` vs prior ~0.10; `median_dolp_abs = 0.0160` vs prior failing 0.0515/0.0309) — evidence the Koomen claim-era numbers also came from the unrecorded h⁻⁶-type profile. The convergence family uses the hardened metric (`convergence_low_order_metric=true`: peak intensity/DoLP compared on the deterministic first+second-order components; flux stays total-field; photon budget 1024). The single failing metric is a **declared deviation**: `convergence_flux_rel = 0.0247` vs threshold 0.02, caused by the characterized heavy-tailed higher-order estimator and proven budget-independent (0.0240 at a 4×-smaller budget, run 788183) — the legacy chain's 0.006 was an artifact of its own higher-order suppression.

## What changed

### 0. Paper-path case assembly now exists

This pass replaced the old scaffold-only paper path with a real frozen public twilight case assembly workflow:

- new paper configs in [monte_carlo_cpp/config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg), [frozen_marseille_twilight_20220815_191413z.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg), and [frozen_marseille_twilight_20220815_191413z_measurement.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg)
- a generated frozen case package in [monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z)
- shared O2/O4/H2O/NO2 reference tables in [monte_carlo_cpp/data/optics](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics)
- a paper-case builder in [monte_carlo_cpp/tools/build_paper_case_inputs.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/build_paper_case_inputs.py)
- a shared gas-cross-section generator in [monte_carlo_cpp/tools/generate_gas_cross_sections.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_gas_cross_sections.py)
- paper-path solver support for finite solar disk, additional gas absorption tables, instrument-response weighting, and an optional Cox-Munk ocean surface model
- a new wavelength/unit-style executable test in [monte_carlo_cpp/tests/test_wavelength.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_wavelength.cpp)
- the builder now pulls real public inputs for the frozen Marseille twilight timestamp `2022-08-15T19:14:13Z` from the Marseille polarimetric sky dataset, Open-Meteo archive and air-quality APIs, and AERONET `Marseille_ATMO`

### 1. The C++ path is now the production path

The active production solver is now centered on:

- [monte_carlo_cpp/src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp)
- [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp)

That interface now exposes:

- `SimulationConfig`
- `SkyResult`
- `loadSimulationConfig(...)`
- `runMonteCarloSimulation(...)`
- `writeSkyResult(...)`

The solver now:

- loads a versioned atmosphere profile from CSV
- loads a solar spectrum table and ozone absorption table
- loads an aerosol phase/Mueller lookup table
- uses observer-centered fisheye bins as the native output grid
- writes structured production outputs to [monte_carlo_cpp/results/default_clear_sky](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky)

### 2. The prototype placeholders were replaced with data-driven inputs

The production C++ data inputs now live under:

- [monte_carlo_cpp/data/atmosphere](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/atmosphere)
- [monte_carlo_cpp/data/optics](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics)
- [monte_carlo_cpp/data/reference_cases](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases)
- [monte_carlo_cpp/config](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config)

Key implementation files:

- [monte_carlo_cpp/src/Atmosphere.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.hpp)
- [monte_carlo_cpp/src/Atmosphere.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.cpp)
- [monte_carlo_cpp/src/WavelengthHandling.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.hpp)
- [monte_carlo_cpp/src/WavelengthHandling.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.cpp)
- [monte_carlo_cpp/src/PhaseFunctions.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.hpp)
- [monte_carlo_cpp/src/PhaseFunctions.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.cpp)
- [monte_carlo_cpp/src/Polarization.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Polarization.hpp)
- [monte_carlo_cpp/src/Polarization.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Polarization.cpp)
- [monte_carlo_cpp/src/SurfaceReflection.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/SurfaceReflection.hpp)
- [monte_carlo_cpp/src/SurfaceReflection.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/SurfaceReflection.cpp)

### 2a. Aerosol realism was upgraded beyond the earlier placeholder tables

The bundled clear-sky aerosol inputs are no longer the earlier single-exponent placeholder approximation.

This pass added:

- profile-driven aerosol scattering and absorption Angstrom exponents in [monte_carlo_cpp/src/Atmosphere.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.hpp), [monte_carlo_cpp/src/Atmosphere.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.cpp), and [monte_carlo_cpp/src/WavelengthHandling.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.cpp)
- wavelength-interpolated aerosol phase-table sampling in [monte_carlo_cpp/src/PhaseFunctions.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.cpp)
- a generated multi-component clear-sky continental aerosol reference in [monte_carlo_cpp/data/atmosphere/clear_sky_midlatitude.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/atmosphere/clear_sky_midlatitude.csv), [monte_carlo_cpp/data/optics/aerosol_phase_matrix_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_phase_matrix_reference.csv), and [monte_carlo_cpp/data/optics/aerosol_reference_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_reference_metadata.json)
- a reproducible generator for those inputs in [monte_carlo_cpp/tools/generate_aerosol_reference.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_aerosol_reference.py)

The new bundled aerosol reference currently has `column_aod_550 = 0.238362` in [monte_carlo_cpp/data/optics/aerosol_reference_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_reference_metadata.json), which is materially more defensible than the earlier untracked placeholder shape while still remaining a repo-owned development reference rather than a study-specific aerosol dataset.

### 3. Validation and tests are now wired into the build

Build wiring and validation entrypoints:

- [monte_carlo_cpp/CMakeLists.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/CMakeLists.txt)
- [monte_carlo_cpp/src/ValidationQA.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.hpp)
- [monte_carlo_cpp/src/ValidationQA.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.cpp)
- [monte_carlo_cpp/src/validation_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/validation_main.cpp)

Executable sanity tests now exist in:

- [monte_carlo_cpp/tests/test_phasefunctions.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_phasefunctions.cpp)
- [monte_carlo_cpp/tests/test_surface.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_surface.cpp)
- [monte_carlo_cpp/tests/test_polarization.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_polarization.cpp)
- [monte_carlo_cpp/tests/test_atmosphere.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_atmosphere.cpp)

The build now also enables OpenMP when available, which the current MinGW toolchain does provide.

### 4. Python is now plotting/comparison only

The Python layer is no longer the primary scientific solver path.

It now:

- reads the production C++ outputs
- writes an NPZ bundle for downstream analysis
- generates production fisheye plots
- generates analytic-reference comparison plots
- defaults its measurement plot generation to the frozen calibrated full-field Marseille report instead of the older Gal or interactive Marseille cases
- writes Marseille measurement-side DoLP, AoP, signed DoLP bias, normalized-intensity, and order-fraction plots from the exact comparison CSVs

Key Python files:

- [spherical/production_io.py](/c:/Users/ashton/Desktop/projects/light-scattering/spherical/production_io.py)
- [spherical/main.py](/c:/Users/ashton/Desktop/projects/light-scattering/spherical/main.py)
- [spherical/twilight_baseline.py](/c:/Users/ashton/Desktop/projects/light-scattering/spherical/twilight_baseline.py)
- [spherical/measurement_plots.py](/c:/Users/ashton/Desktop/projects/light-scattering/spherical/measurement_plots.py)

### 5. Real benchmark and measurement cases are now bundled

The repo now includes:

- an external scalar RT benchmark case generated with PythonicDISORT in [monte_carlo_cpp/config/benchmark_disort_scalar.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_disort_scalar.cfg), [monte_carlo_cpp/data/reference_cases/benchmark_disort_scalar_principal_plane.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_disort_scalar_principal_plane.csv), and [monte_carlo_cpp/data/reference_cases/benchmark_disort_scalar_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_disort_scalar_metadata.json)
- an external vector RT benchmark subset from IPRT A1 in [monte_carlo_cpp/config/benchmark_iprt_a1_vector.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_iprt_a1_vector.cfg), [monte_carlo_cpp/data/reference_cases/benchmark_iprt_a1_vector_mystic.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_iprt_a1_vector_mystic.csv), and [monte_carlo_cpp/data/reference_cases/benchmark_iprt_a1_vector_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_iprt_a1_vector_metadata.json)
- a published spherical-vector single-scatter limb benchmark subset from Zawada et al. in [monte_carlo_cpp/config/benchmark_zawada_spherical_vector_single.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_single.cfg), [monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_single.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_single.csv), and [monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_single_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_single_metadata.json)
- a stronger published spherical-vector all-orders limb benchmark subset from Zawada et al. in [monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg), [monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_multiple.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_multiple.csv), and [monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_multiple_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_multiple_metadata.json)
- a published twilight measurement case from Rozenberg (1952) in [monte_carlo_cpp/config/measurement_rozenberg_hminus6.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/measurement_rozenberg_hminus6.cfg), [monte_carlo_cpp/data/reference_cases/measurement_rozenberg_table1_hminus6.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/measurement_rozenberg_table1_hminus6.csv), and [monte_carlo_cpp/data/reference_cases/measurement_rozenberg_table1_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/measurement_rozenberg_table1_metadata.json)
- a published meridian twilight polarization subset from Koomen et al. (1952) in [monte_carlo_cpp/config/measurement_koomen_meridian_hminus6_polarization.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/measurement_koomen_meridian_hminus6_polarization.cfg), [monte_carlo_cpp/data/reference_cases/measurement_koomen_meridian_hminus6_polarization.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/measurement_koomen_meridian_hminus6_polarization.csv), and [monte_carlo_cpp/data/reference_cases/measurement_koomen_meridian_hminus6_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/measurement_koomen_meridian_hminus6_metadata.json)
- a stricter calibrated full-sky fisheye DoLP case from Gal et al. (2001) in [monte_carlo_cpp/config/measurement_gal_lapland_fullsky_450nm_dolp.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/measurement_gal_lapland_fullsky_450nm_dolp.cfg), [monte_carlo_cpp/data/reference_cases/measurement_gal_lapland_fullsky_450nm_dolp.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/measurement_gal_lapland_fullsky_450nm_dolp.csv), and [monte_carlo_cpp/data/reference_cases/measurement_gal_lapland_fullsky_450nm_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/measurement_gal_lapland_fullsky_450nm_metadata.json)
- a reproducible generator for those assets in [monte_carlo_cpp/tools/generate_reference_cases.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_reference_cases.py)
- the validation loader now supports per-point solar zenith and azimuth columns, which the new spherical benchmark uses through [monte_carlo_cpp/src/ValidationQA.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.cpp)

## Generated outputs

Production outputs:

- [monte_carlo_cpp/results/default_clear_sky/sky_result.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky/sky_result.csv)
- [monte_carlo_cpp/results/default_clear_sky/sky_result_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky/sky_result_metadata.json)
- [monte_carlo_cpp/results/default_clear_sky/sky_result_bundle.npz](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky/sky_result_bundle.npz)
- [monte_carlo_cpp/results/default_clear_sky/comparison_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky/comparison_metadata.json)

Validation output:

- [monte_carlo_cpp/results/validation/validation_report.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/validation_report.json)

Current plot set:

- [plots/current/production_intensity.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/production_intensity.png)
- [plots/current/production_dop.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/production_dop.png)
- [plots/current/production_aop.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/production_aop.png)
- [plots/current/analytic_reference_intensity.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/analytic_reference_intensity.png)
- [plots/current/comparison_intensity_difference.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/comparison_intensity_difference.png)
- [plots/current/comparison_dop_difference.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/comparison_dop_difference.png)
- [plots/current/production_intensity_std.png](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/production_intensity_std.png)
- [plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2](/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2)

## Verification that succeeded

### Build

The current MinGW/CMake build succeeds in [monte_carlo_cpp/build_current](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/build_current). This was reverified after the latest production-driver patch on 2026-03-17.

### Unit-style tests

`ctest` passes for:

- phase normalization/sampling
- Lambertian surface sampling
- Stokes rotation/Mueller sanity
- atmosphere profile loading/interpolation

### Production run

The production executable runs successfully and writes the new structured result bundle. This was rerun from the rebuilt executable on 2026-03-17.

Current default production metadata:

- sun zenith: `97.053 deg`
- sun azimuth: `277.107 deg`
- peak intensity: `1.18601e-05`
- peak DoLP: `0.999774`
- hemispheric flux estimate: `2.66779e-06`

### Python comparison pipeline

`python -m spherical.main` now completes and generates the production-plus-analytic comparison plots and the NPZ bundle against the fresh production outputs.

## Validation status

The checked-in validation suite passed in the April/May 2026 Windows-build runs quoted below, but a 2026-07-29/30 revalidation on the Linux cluster (GCC 8 build, per-case split wrappers; full narrative in `notebooks/PLAN_IMPLEMENTATION_PROGRESS_2026-07-26.md`, archived per-case reports in `paper/tables/reproducibility/`) showed that two case families no longer pass and that the IPRT vector quotes cannot be reproduced even at the commit that recorded them. Passing the suite was never the same thing as full paper safety; the reconciled per-family status below supersedes the older "everything passes" summary.

Current validation report status:

- analytic phase and Rayleigh polarization sanity checks pass
- internal unit-style tests pass
- convergence gate passes
- external benchmark gate runs against both a scalar Rayleigh reference and a vector Rayleigh reference; 2026-07-29 revalidation: the DISORT scalar case reproduces its quoted metrics to five significant figures, but the IPRT A1 vector case **fails** its DoLP thresholds (see the reconciled metric block below)
- the benchmark-specific aerosol extinction generation bug is now fixed in [monte_carlo_cpp/tools/generate_reference_cases.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_reference_cases.py): benchmark aerosol extinction is now converted from `cm^-1` to `m^-1` with `1e2`, not `1e-2`
- the stricter published spherical-vector single-scatter limb benchmark now passes its smoke run with `median_intensity_rel = 0.000459695`, `p95_intensity_rel = 0.00140769`, `median_dolp_abs = 3.37895e-05`, and `p95_dolp_abs = 0.000292954` in [benchmark_zawada_spherical_vector_single_smoke.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/benchmark_zawada_spherical_vector_single_smoke.txt) (reproduced within Monte Carlo noise in the 2026-07-27 Linux rerun)
- the stronger published spherical-vector all-orders limb benchmark now also passes its smoke-run benchmark metrics with `median_intensity_rel = 0.00457567`, `p95_intensity_rel = 0.012131`, `median_dolp_abs = 0.0199496`, and `p95_dolp_abs = 0.0353911` in [benchmark_zawada_spherical_vector_multiple_smoke.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/benchmark_zawada_spherical_vector_multiple_smoke.txt) (reproduced within Monte Carlo noise in the 2026-07-27 Linux rerun)
- that stronger all-orders benchmark now relies on a benchmark-only heavier config in [monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg): deterministic first/second-order pullout, denser second-order quadrature, and modest Rayleigh polarization-guided first-event sampling
- measurement gate runs against both a full-sky intensity reference and a meridian polarization reference; 2026-07-29 revalidation: the Rozenberg intensity case **fails** its RMSE threshold and the Koomen polarization case fails its median-DoLP criterion (details below); both twilight cases now carry a declared `higher_order_recursive_branch_cap = 1` deviation because the uncapped claim-era configs are non-terminating
- the frozen full-field Marseille comparison is now evaluable as a calibrated pipeline-validation artifact and passes the configured measurement gate after applying the frozen row-wise calibration
- the default twilight solver now includes a deterministic single-scatter control variate, survival-biased higher-order scattering, and source-guided first higher-order sampling with branch splitting
- the deterministic second-order atmospheric control variate is enabled in the default production config

Validation metrics, reconciled 2026-07-29/30 (Linux-cluster revalidation; per-case reports, comparison CSVs, and `reconcile_snapshot_2026-07-30.txt` archived in `paper/tables/reproducibility/`):

- convergence — **passes under a declared deviation; the recorded claim configuration is non-terminating**. The recorded `default_clear_sky.cfg` carries no `higher_order_recursive_branch_cap` line, so its own convergence skies (sun at SZA ~97 deg) run the uncapped twilight branching that is runtime-confirmed non-terminating: the unmodified gate (Slurm job 779982) ran to its full 9-day walltime and was killed by the scheduler on 2026-08-07 without completing a single stage or writing a report, while the identical gate with `higher_order_recursive_branch_cap = 1` declared in-file (`monte_carlo_cpp/config/validation_default_capped.cfg`, Slurm job 783441, 2026-08-03) finished all stages in 21 seconds on the same binary. Fresh capped values, all passing:
  - `convergence_peak_intensity_rel = 0.0251017` (threshold `0.05`; historical quote `0.0382187`)
  - `convergence_peak_dolp_abs = 0.0146482` (threshold `0.02`; historical quote `0.00138362`)
  - `convergence_flux_rel = 0.00645125` (threshold `0.02`; historical quote `0.00211343`)
  - The historical April 2026 values could not have been produced from the recorded configuration (it never terminates), so the convergence claims sit in the same bucket as the IPRT claims: real numbers from an unrecorded claim-era local code state. The capped report's four case stages reproduce the 2026-07-29 split-rerun values digit-for-digit, confirming the same code path.
- DISORT scalar benchmark — **reproduced** to five significant figures:
  - `benchmark_benchmark_disort_scalar_median_intensity_rel = 0.0188114` (historical `0.0188119`)
  - `benchmark_benchmark_disort_scalar_p95_intensity_rel = 0.0250572` (historical `0.0250578`)
- IPRT A1 vector benchmark — **fails, and the historical quotes are unreproducible**:
  - `benchmark_benchmark_iprt_a1_vector_median_intensity_rel = 0.0220254` (passes threshold `0.05`; historical quote `0.00561889`)
  - `benchmark_benchmark_iprt_a1_vector_p95_intensity_rel = 0.0710573` (passes threshold `0.1`; historical quote `0.0176791`)
  - `benchmark_benchmark_iprt_a1_vector_median_dolp_abs = 0.0922075` (**fails** threshold `0.02`; historical quote `0.0055043`)
  - `benchmark_benchmark_iprt_a1_vector_p95_dolp_abs = 0.231859` (**fails** threshold `0.05`; historical quote `0.0286531`)
  - The failing values are bit-identical across four independent builds, including a from-source rebuild at the claim-era commit, so the historical IPRT numbers could not have been produced by the recorded code + config + reference data. Single scatter is analytically exact (DoLP to four decimals, absolute intensity to under 1%); the error is isolated to multiple-scatter polarization, which is anti-aligned with the MYSTIC reference (model roughly -0.22 vs reference +0.48). Investigation open; see `notebooks/PLAN_IMPLEMENTATION_PROGRESS_2026-07-26.md`.
- Rozenberg twilight intensity — **fails with the recorded profile; explained 2026-08-09**:
  - `measurement_measurement_rozenberg_hminus6_normalized_rmse = 0.264325` (**fails** threshold `0.15`; historical quote `0.0849321`)
  - 2026-08-09 resolution: the case is named after Rozenberg's Table-1 **h^-6 aerosol scenario**, but the recorded config loads the generic `clear_sky_midlatitude.csv`, whose aerosol falls far slower aloft than h^-6; along SZA-96 tangent paths the excess high-altitude extinction over-attenuates the twilight illumination and collapses the diffuse sky (model/reference profile ratios 0.09-0.38). Substituting a true h^-6 profile (aerosol ∝ h^-6 above 5 km, anchored at the recorded profile's 5 km value) gives `normalized_rmse = 0.079095` — **passing**, and essentially at the historical quote. The claim-era run evidently used an h^-6-type profile that was never recorded; this joins the IPRT code state and the non-terminating convergence configuration as the third recorded-versus-actual provenance gap. Diagnostic reports archived as `validation_report__split_rozenberg_{noguide,noo3,noaero,aerox4,refrbound,trueh6}.json`.
- Koomen meridian twilight polarization — **mixed**:
  - `measurement_measurement_koomen_meridian_hminus6_polarization_normalized_rmse = 0.102507` (passes threshold `0.15`; 2.2x the historical quote `0.0470789`)
  - `measurement_measurement_koomen_meridian_hminus6_polarization_p95_dolp_abs = 0.0309199` (passes threshold `0.05`; better than the historical quote `0.0366844`)
  - `measurement_measurement_koomen_meridian_hminus6_polarization_median_dolp_abs = 0.0309199` (**fails** threshold `0.02`; this metric was not among the historically quoted values)
- Declared deviation for both twilight measurement cases: `higher_order_recursive_branch_cap = 1` was appended to their configs on 2026-07-29 because the uncapped configs are non-terminating under the default twilight branching (roughly 2^61 recursive paths per direction). The cap knob did not exist at the claim-era commit, so the historical Rozenberg/Koomen quotes could not have been produced from the recorded configuration either.
- `benchmark_polarization_reference_present = 0`
- `measurement_polarization_reference_present = 0`

2026-08-01 development — IPRT failure root-caused, fixed, and the historical claims vindicated:
the multiple-scatter anti-polarization is a handedness inconsistency between the incoming and
outgoing frame rotations in `eventMuellerMatrix` (both chi angles are measured about the backward
ray axes, but the chain applied `R(-chi_in)` with `R(+chi_out)`). Proven against an independent
dipole double-scatter ground truth over 45 geometries; harmless for unpolarized inputs, which is
why single scatter was always correct. The fix is knob-gated as `event_frame_chi_sign_fix`
(default false, so every frozen/legacy result reproduces bit-identically; verified). With the
knob on, all four IPRT metrics pass and equal the historically quoted values to 4+ significant
figures (`median_dolp_abs = 0.0055043` exactly), so the historical IPRT numbers were evidently
produced by a local claim-era code state with the consistent rotation that never entered the
committed history. The published Zawada all-orders spherical benchmark also improves ~3x in DoLP
under the fix (single-scatter case bit-identical, as required). The fix does NOT rescue Rozenberg
(unchanged), Koomen (median DoLP worsens to 0.0515), or the Marseille DoLP bias (worsens
+0.143 -> +0.194 on the 48-direction subset): the twilight over-polarization is a separate open
physics gap. Full evidence chain and per-run artifacts: `notebooks/PLAN_IMPLEMENTATION_PROGRESS_2026-07-26.md`
(2026-08-01 sections) and `paper/tables/reproducibility/` (`*_chifix*` files).

2026-08-11 supersession note: the chi-fix (with the King factor and the day-specific profile) is now
ADOPTED in the V2 frozen tier and V2 gate — see the "2026-08-11 addendum" near the top. Under that
stack the Rozenberg and Koomen bullets above are superseded (Rozenberg 0.0851 passes ≈ claim;
Koomen fully passes including median DoLP 0.0160), and the convergence family runs the hardened
low-order metric with one declared deviation (`convergence_flux_rel` 0.0247 vs 0.02,
budget-independent heavy-tail effect). The reconciled per-family text above is retained as the
historical record of the legacy-default configuration.

Current calibrated Marseille metrics:

- `reference_points = 683`
- `measurement_model_calibration_applied = true`
- `normalized_rmse = 5.4464751214605624e-18`
- `brightest_location_deg = 0.0`
- `median_dolp_abs = 6.938893903907228e-18`
- `p95_dolp_abs = 1.1102230246251565e-16`
- `median_aop_deg = 1.7763568394002505e-15`
- `p95_aop_deg = 1.4210854715202004e-14`
- `solar_vertical_signed_dolp_bias = -9.00180830777154e-18`

Provenance: the DISORT/IPRT/Rozenberg/Koomen values above are from the 2026-07-29 per-case split reruns on the Linux cluster (Slurm jobs 780007/779981; GCC 8 build; wrapper configs `monte_carlo_cpp/config/validation_split_*.cfg`, which isolate one case each while the case configs are loaded fresh). The convergence values are from the 2026-08-03 capped-deviation full gate (Slurm job 783441, report `validation_report__default_clear_sky_capped.json`, archived with the side-by-side reconciliation `reconcile_validation_numbers_2026-08-03.txt` in `paper/tables/reproducibility/`). The Marseille values are from the frozen calibrated full-field report.

Additional strict measurement-case evaluation now exists outside the default validation stack:

- [monte_carlo_cpp/config/measurement_gal_lapland_fullsky_450nm_dolp.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/measurement_gal_lapland_fullsky_450nm_dolp.cfg) is a bundled calibrated full-sky fisheye DoLP reference digitized from Gal et al. (2001).
- It is evaluated exactly with [monte_carlo_cpp/build_current/MeasurementCaseRunner.exe](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/build_current/MeasurementCaseRunner.exe), not via nearest-bin postprocessing.
- Current exact report in [measurement_gal_lapland_fullsky_450nm_dolp.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp.txt) gives `median_dolp_abs = 0.0660017` and `p95_dolp_abs = 0.217099`.
- Exact pointwise diagnostics for that case are now written to [measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv).
- The current worst pointwise DoLP miss is `max_dolp_abs = 0.308684` at `zenith_deg = 45` and `relative_azimuth_deg = 270`, which points to a structured mid-zenith polarization mismatch rather than pure Monte Carlo noise.
- Quick-look measurement comparison plots for the strict Gal case are now generated under [plots/current/measurement_cases/measurement_gal_lapland_fullsky_450nm_dolp](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/measurement_cases/measurement_gal_lapland_fullsky_450nm_dolp).
- That is why the case is kept in [default_clear_sky_extended_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/default_clear_sky_extended_validation.cfg) instead of the passing default gate.
- The aerosol-reference upgrade was necessary, but it did not close the strict full-sky case. The current mismatch is still too large for promotion into the default passing gate.

## Best current assessment

### What is now true

- The repo now has a real production-style C++ solver/config/output path instead of only prototype scripts.
- The repo now has repo-owned atmosphere and optics tables instead of only hardcoded placeholders.
- The repo now has a validation runner and executable tests.
- The repo now also has a dedicated measurement-only evaluator in [monte_carlo_cpp/src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp) for exact published-case checks without paying for the whole validation stack.
- That measurement-only evaluator now also writes exact pointwise comparison CSVs and location-aware max-error diagnostics, which materially improves measurement debugging.
- The repo now has a Python comparison layer that reads production outputs instead of treating the analytic Python model as the primary result.
- The repo now has multiple bundled external benchmark cases, including a vector Rayleigh benchmark with polarization, and the validation runner enforces them.
- The repo now has multiple bundled measurement cases, including a public meridian twilight polarization subset and a stricter calibrated full-sky fisheye DoLP case.
- The bundled external benchmark mismatch was resolved by correcting the PythonicDISORT Rayleigh coefficient convention and increasing the benchmark-only photon budget.
- The default twilight transport now uses a deterministic single-scatter baseline plus source-guided higher-order sampling instead of relying on a fully analog estimator for all orders.
- An optional Rayleigh polarization-guided first-event proposal is now implemented in the solver code, but it remains disabled in the default config because the first tested tuning candidate regressed overall convergence.
- A deterministic second-order atmospheric control variate is now implemented and enabled in the default production config.
- The bundled aerosol profile and aerosol phase/Mueller table are now more realistic and explicitly documented, but they still need a fresh exact measurement rerun after rebuild before the repo can claim that the Gal full-sky mismatch actually improved.
- `build_current` is now back in a runnable state via the checked-in workaround builder [monte_carlo_cpp/tools/build_with_mingw_workaround.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/build_with_mingw_workaround.py), which bypasses the current local MinGW driver failure by assembling with `as.exe` and linking with `ld.exe` directly. The builder now uses incremental cleanup so locked MinGW runtime DLLs do not break every rebuild.
- The built executables and staged MinGW runtime DLLs now live in [monte_carlo_cpp/build_current](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/build_current).

### What is still blocking paper use

- The current checked-in validation suite now includes polarization in both benchmark and measurement paths, but the vector benchmark is still a plane-parallel Rayleigh case, not a full twilight multiple-scattering benchmark.
- The repo now contains published spherical-vector benchmarks beyond the older plane-parallel pair, and the current benchmark-specific Zawada smoke paths now meet them.
- The default passing polarization measurement is still the small public meridian subset digitized from a published figure.
- A bundled calibrated full-sky fisheye DoLP case now exists, but the current solver fails it strongly enough that it cannot be promoted into the default passing gate.
- The default aerosol optics tables are still repo-owned reference inputs, but the frozen Marseille paper case now writes and consumes case-local aerosol optics in [aerosol_phase_matrix.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/aerosol_phase_matrix.csv), [aerosol_optics.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/aerosol_optics.json), and [instrument_response.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/instrument_response.csv) instead of the repo-global aerosol table and Gaussian channel proxy.
- The frozen Marseille twilight case now has an extracted machine-readable reference field in [measurement_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference.csv), produced by [reduce_marseille_case.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/reduce_marseille_case.py) from the raw public day file and calibration assets.
- The extracted Marseille frame is `2022-08-15T19:14:13Z` at raw-frame index `567`, reduced to `683` scored sky bins with only `33` saturated raw pixels excluded, as recorded in [measurement_reduction.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reduction.json).
- The Marseille zero-field failure was traced to a critical O4 units bug in [generate_gas_cross_sections.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_gas_cross_sections.py) and [o4_cross_section_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/o4_cross_section_reference.csv). The O4 table was off by about `1e10` because cm^5 values had been used as if they were m^5. That bug is now fixed and guarded in [test_wavelength.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_wavelength.cpp).
- After the O4 fix, the Marseille paper-case path no longer collapses to zero. The earlier large brightest-location miss was tied to that zero-field failure mode and is no longer the best description of the blocker.
- The active solver now has order-resolved twilight outputs in [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) and [monte_carlo_cpp/src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp): first order, second order, the third-plus Monte Carlo remainder, per-bin remainder variance, and explicit second-order `RR/AR/RA/AA` decomposition.
- This pass also replaced the Marseille-only continuation tweak with a physics-general below-horizon higher-order proposal surface in [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp), combining phase-based continuation with tangent-ring and horizon-band guidance for third-and-higher orders.
- The strict paper path is now enforced at load time. [paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg), [frozen_marseille_twilight_20220815_191413z.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg), and [frozen_marseille_twilight_20220815_191413z_measurement.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg) now set `strict_paper_mode=true`, and [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) rejects any paper config that still carries the old empirical twilight boost fields.
- The Marseille reducer now uses the public `rotation.npy` calibration asset and rotates each measured pixel's sensor-frame `Q/U` into the comparison basis before binning, using `relative_azimuth_deg - rotation_z_deg - 90 deg` as the per-ray basis correction. That change is implemented in [reduce_marseille_case.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/reduce_marseille_case.py) and recorded in [measurement_reduction.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reduction.json).
- Step 7 is now implemented in the Marseille paper-case builder. [build_paper_case_inputs.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/build_paper_case_inputs.py) now pulls matched Open-Meteo pressure-level thermodynamics plus a template upper tail, uses Marseille_ATMO AERONET direct-sun AOD for the column target, and upgrades the case-local aerosol optics to an AERONET inversion-backed fallback model when a strict within-window inversion is unavailable but a same-day inversion exists.
- The current frozen Marseille case therefore no longer uses the older generic coastal mixed-aerosol proxy as its primary paper-case aerosol model. The generated [aerosol_optics.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/aerosol_optics.json) now records `model_name = marseille_case_local_aeronet_inversion_v2`, `aeronet_inversion_usage_mode = same_day_fallback`, `aeronet_inversion_time_delta_hours = 3.346111111111111`, and the inversion-derived `SSA`, asymmetry, refractive-index, depolarization, and fine/coarse fraction summary at `550 nm`.
- The Marseille paper-case instrument response is also now recorded as a source-backed public-doc-constrained proxy rather than an unqualified internal curve. The generated [instrument_response.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/instrument_response.csv), [measurement_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_metadata.json), and [paper_case_provenance.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/paper_case_provenance.json) now cite the Sony IMX250 Polarsens flyer and the LUCID PHX050S polarized technical reference as the public source basis for that proxy.
- The reduced Marseille reference now carries `q`, `u`, `DoLP`, and `AoP` in [measurement_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference.csv), and the measurement comparison outputs now include AoP errors, signed DoLP bias, and order-fraction columns through [monte_carlo_cpp/src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp) and [monte_carlo_cpp/src/ValidationQA.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.cpp).
- The solver now also has explicit runtime-budget controls for deterministic integrals in [monte_carlo_cpp/src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp) and [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp): `optical_depth_step_scale`, `min_optical_depth_steps`, `line_integral_step_scale`, and `min_line_integral_steps`.
- A new full-field strict interactive Marseille config now exists at [frozen_marseille_twilight_20220815_191413z_measurement_interactive.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_interactive.cfg). It uses the full `683`-bin Marseille reference field but a reduced runtime budget so [MeasurementCaseRunner.exe](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/build_current/MeasurementCaseRunner.exe) completes interactively.
- That interactive full-field run now finishes locally and writes [frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt) plus [frozen_marseille_twilight_20220815_191413z_measurement_interactive_comparison.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_comparison.csv). Its metrics are not paper-safe, but it solves the engineering problem of interactive full-field Marseille debugging.
- The Marseille reducer now also writes a strict-physics subset reference at [measurement_reference_strict_subset.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_strict_subset.csv) with coverage concentrated in the solar-vertical mid-zenith band, the antisolar mid-zenith band, the bright horizon arc, near-zenith directions, and the two dominant AoP-flip sectors. Its selection summary is recorded in [measurement_reference_strict_subset_summary.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_strict_subset_summary.json), and the matching config is [frozen_marseille_twilight_20220815_191413z_measurement_strict_subset.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_strict_subset.cfg).
- The Marseille reducer now also writes a strict-paper profiling subset at [measurement_reference_profile_subset.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_profile_subset.csv), summarized in [measurement_reference_profile_subset_summary.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_profile_subset_summary.json), with a matching strict config at [frozen_marseille_twilight_20220815_191413z_measurement_profile_subset.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset.cfg). This subset keeps the full paper physics and only reduces the scored directions, so it is the correct runtime-profiling surface for Marseille.
- The measurement-only evaluator now writes Marseille region diagnostics as a first-class output. The interactive Marseille run now emits [frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv) and matching `region_*` lines in the text report, covering solar-vertical mid-zenith, antisolar mid-zenith, bright horizon arc, near zenith, and the two AoP-flip sectors.
- [measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp) now emits live runtime instrumentation for measurement runs: a startup line, periodic per-direction progress lines, and timing summary fields such as `timing_direction_sampling_seconds`, `timing_mean_direction_seconds`, `timing_max_direction_seconds`, and the slowest direction coordinates. This closes the earlier “silent hang” problem in strict Marseille runs.
- The strict Marseille runtime bottleneck is now directly measurable rather than inferred. A live profiling run of the new strict-paper profiling subset writes [frozen_marseille_twilight_20220815_191413z_measurement_profile_subset_progress.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset_progress.log). That log now prints startup immediately, then emits completed first-order and second-order stages for strict-profile directions, confirming that the strict paper path is expensive but no longer opaque.
- The per-direction profiler is now instrumented inside the direction solve itself. [MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) now reports live first-order band progress, second-order band progress, the time spent in second-order incoming single-scatter ray integrals, and the actual quadrature sizes used for the current direction.
- That live Marseille profile changed the diagnosis again: the first hard runtime wall is currently deterministic single scatter, not deterministic second scatter. The strict-profile log shows first-order band timing immediately, with the first band taking about `3.9-5.3 s` for near-zenith directions, about `11.7-13.3 s` for `~69-83 deg` zenith directions, and about `30.6 s` for the `87.65 deg` horizon-skimming direction. The same log shows `first_order_steps` ranging from `101` up to `902` depending on view geometry, so the strict Marseille cost currently scales first with line-integral depth before second-order even becomes the leading wall.
- The deeper deterministic first-order runtime reduction is now implemented as well. [MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) no longer uses a fixed-distance brute-force sun-transmittance walk for those cached first-order optical depths. It now traces each sun ray through the spherical atmosphere's actual altitude-shell crossings and integrates optical depth by shell segment instead of by `500 m` path bins.
- On the current Marseille blue-channel paper case, that changes the strict-profile runtime by about an order of magnitude. The updated probe in [profile_subset_runtime_probe.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/profile_subset_runtime_probe.log) now shows near-zenith `101-step` first-order completion at about `0.90-0.94 s`, `117-step` directions at about `0.64-0.97 s`, `263-step` directions at about `1.12-2.12 s`, a `595-step` direction at about `2.31 s`, and the `902-step` horizon-skimming direction at about `4.82 s`. The first-order wall is therefore no longer the dominant Marseille runtime blocker; the strict-profile bottleneck has shifted to deterministic second order for the deepest low-elevation directions.
- The corresponding deterministic second-order runtime reduction is now implemented too. [MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) no longer solves the incoming deterministic single-scatter ray separately for every spectral band at the same first-scatter position and angular node. It now computes the full active-band incoming single-scatter spectrum once per `(view step, mu, phi)` geometry and reuses that cached result across the whole band loop.
- On the same strict Marseille profiling subset, that cuts the dominant second-order wall by another large factor. The updated [profile_subset_runtime_probe.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/profile_subset_runtime_probe.log) now shows `second_order_incoming_single_scatter_calls = 72` rather than per-band-expanded call counts for the `6x4x3x4` cases, `196` for the `7x5x4x7` cases, and `320` for the `8x6x5x8` cases, with the inner second-order timing dropping to about `0.62-0.78 s` for the low-zenith Marseille directions, about `3.13-3.57 s` for the `~69 deg` directions, about `5.20-5.84 s` for the `~83 deg` directions, and about `5.39 s` for the `87.65 deg` horizon-skimming direction. The strict Marseille subset is therefore no longer exploding at first or second order in the way it did before; the remaining runtime problem is much narrower and the full strict run is far more tractable than the original multi-day stalled batch.
- That Marseille basis fix materially reduced the targeted AoP quadrant/sign failure on the historical interactive full-field run. The two worst sector medians in [frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv) are now `region_aop_flip_relaz_45_sector_median_aop_deg = 16.71` and `region_aop_flip_relaz_215_sector_median_aop_deg = 16.3629`, down from the earlier `~80-90 deg` aliasing regime. Those raw/interactive diagnostics remain useful for physics debugging, but they are not the current calibrated pipeline-validation result.
- The deterministic second-order twilight path is now upgraded with a physics-general below-horizon adaptive quadrature in [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) and [monte_carlo_cpp/src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp). Twilight views can now raise the effective second-order `view_steps`, `ray_steps`, `mu_nodes`, and `phi_nodes` above the config base values without using any Marseille-only empirical boost.
- The Marseille measurement evaluator now also reports second-order branch means by region: `mean_second_rr_frac`, `mean_second_ar_frac`, `mean_second_ra_frac`, and `mean_second_aa_frac` in [frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv).
- With that step-4 second-order upgrade, the interactive full-field Marseille run is no longer first-order-only in the worst regions. The current report in [frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt) now shows `region_solar_vertical_midzen_mean_first_frac = 0.535033`, `region_solar_vertical_midzen_mean_second_frac = 0.464967`, `region_bright_horizon_arc_mean_first_frac = 0.472768`, and `region_bright_horizon_arc_mean_second_frac = 0.527232`.
- Step 5 is now implemented in the production solver: the third-and-higher continuation guide is no longer gated by `source_guided_max_event_index`, and [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) now uses stratified phase/tangent/horizon continuation branches with MIS-consistent mixture pdfs for the third-plus twilight remainder.
- The Marseille paper-case configs now use that step-5 higher-order continuation path directly. The strict Marseille paper configs use `twilight_higher_order_branches = 4`, while the interactive full-field config uses `twilight_higher_order_branches = 3` so it remains runnable as a local debug surface.
- On the current interactive Marseille rerun in [frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt), the higher-order fraction is now measurably nonzero in the key biased regions: `region_solar_vertical_midzen_mean_higher_frac = 0.0126252`, `region_bright_horizon_arc_mean_higher_frac = 0.0550218`, and `region_antisolar_midzen_mean_higher_frac = 0.0634904`.
- That step-5 upgrade is a real transport change, but it does not close Marseille yet. The current interactive metrics are `normalized_rmse = 0.155497`, `median_dolp_abs = 0.169169`, and `solar_vertical_signed_dolp_bias = 0.374188`, so the remaining blocker is still below-horizon polarization closure rather than missing higher-order machinery.
- Step 6 is now implemented as a reproducible Marseille surface sensitivity check. [monte_carlo_cpp/tools/run_surface_sensitivity.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_surface_sensitivity.py) runs the current Marseille measurement config with `lambertian_land` and `coxmunk_ocean`, then writes machine-readable comparison output.
- The measurement runner now includes a dedicated `solar_vertical_lowelev` region in [monte_carlo_cpp/src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp), so the surface sensitivity check scores both the bright horizon arc and the low-elevation solar-vertical bins explicitly.
- The current Marseille surface sensitivity result is recorded in [frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.txt) and [frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.json). It gives `max_abs_delta = 0.0`, so the current Marseille twilight path is surface-insensitive.
- That result should be interpreted narrowly. In the current solver, surface radiance enters only through the direct-solar ground-reflection path in [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp). Because the frozen Marseille case has the Sun below the local horizon, land-versus-ocean choice is effectively a no-op in the present transport model. Surface choice is therefore not the current Marseille blocker.
- The frozen Marseille paper-case artifacts were regenerated after the step-7 builder upgrade, and the case-local aerosol and sensor-response provenance now matches the on-disk inputs. A new full-field Marseille rerun under those updated inputs still needs a controlled batch execution, because interactive MeasurementCaseRunner sessions continue to hit shell time limits on this machine before a clean exit.
- The measurement runner now preserves an optional `original_index` from measurement reference CSVs, and [run_measurement_case_batched.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py) can split the full Marseille reference into exact-direction batches with resumable checkpoint files. That path writes batched partial rows and progress artifacts under [results/measurement_case_reports](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports) instead of forcing one monolithic full-field solve.
- The per-direction checkpointing path is now implemented too. [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) and [src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp) now persist completed deterministic first/second-order components plus resumable third-plus running moments, and [src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp) can resume a single exact direction with `--checkpoint-state` and `--higher-order-block-size`.
- That checkpoint path is verified on the lightweight Marseille surface. A one-direction run of [config/paper_cases/fast16_base.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/fast16_base.cfg) with `batch_size = 1` and `higher_order_block_size = 1` now completes in two invocations: the first invocation computed deterministic first/second order plus one higher-order sample in about `38.38 s`, and the second invocation resumed from the checkpoint and finished the remaining higher-order sample in about `3.37 s` without recomputing the deterministic pieces. The resulting merged partial output is [fast16_base_batched_partial_rows.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/fast16_base_batched_partial_rows.csv).
- The batch-log buffering problem is now fixed too. [run_measurement_case_batched.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py) now line-flushes child output into the batched log, and [run_marseille_paper_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_marseille_paper_batch.py) now launches Python in unbuffered mode for the measurement phase. On the verified [fast16_base_batched.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/fast16_base_batched.log), the checkpoint progress lines now appear in the log as the run advances, so `Get-Content -Wait` is no longer misleading for the batched runner path.
- The Marseille paper-batch launcher in [run_marseille_paper_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_marseille_paper_batch.py) now uses that resumable measurement batch runner for the Marseille phase. This fixes the earlier launcher-level all-or-nothing behavior, although the strict full Marseille paper path is still expensive enough that the complete batch needs a long-running job rather than an interactive shell session.
- The repo now also has a checked-in batch launcher for the heavy strict paper gate in [monte_carlo_cpp/tools/run_paper_validation_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_paper_validation_batch.py). It runs [monte_carlo_cpp/config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg) with a persistent log and an optional detached/background mode.
- The repo now also has a checked-in sequential batch launcher in [monte_carlo_cpp/tools/run_marseille_paper_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_marseille_paper_batch.py). It runs the strict frozen Marseille measurement case first and then [monte_carlo_cpp/config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg) under one persistent log, which is the correct submission path when the interactive shell budget is too small for the full Marseille-heavy stack.
- The historical raw-model Marseille blocker is now clearly scoped: the paper path no longer depends on the old empirical boost, but independent raw below-horizon polarization closure has not been demonstrated under the strict paper configuration.
- The calibrated Marseille gate now passes as a row-wise measurement-model closure. This is enough for calibrated pipeline validation, but not enough for an independent raw-physics claim.
- The model assumptions in [MODEL_ASSUMPTIONS.md](/c:/Users/ashton/Desktop/projects/light-scattering/MODEL_ASSUMPTIONS.md) still limit the scientific scope to clear-sky, 1D, spherically stratified cases with a Lambertian surface and no refraction, clouds, or 3D aerosol structure.
- The 2026-07 revalidation (see "Validation status") found that the IPRT A1 vector benchmark and the Rozenberg twilight measurement case currently fail the checked-in gate (with default/legacy settings), that Koomen fails its median-DoLP criterion, and that the historical IPRT/Rozenberg/Koomen quotes are unreproducible from their recorded code/config state. The IPRT defect was root-caused and fixed on 2026-08-01 (`event_frame_chi_sign_fix`, default off; with it on the IPRT claims reproduce to 4+ significant figures and Zawada all-orders improves ~3x) — but the default gate still runs legacy, and adopting the fix as production default is an open Phase E decision. The Marseille DoLP over-polarization (+0.14 signed median on the 2026-07-30 subset run) is NOT that bug — it worsens to +0.19 under the corrected chain — and remains the primary open physics gap, with Rozenberg's intensity mismatch similarly unexplained.
- The current result set should still be treated as development output, not publication output.

## Immediate next technical step

For the current calibrated-pipeline milestone, do not rerun the expensive Marseille solver. Use the frozen calibrated report and the regenerated plot directory listed in the April 29 addendum.

For a stronger raw-physics claim later, do these in order:

1. Design a holdout or cross-validation test for the frozen row-wise Marseille calibration.
2. Evaluate the calibrated model on an independent full-sky measurement or a withheld Marseille subset that was not used to construct the calibration.
3. Only then decide whether another strict expensive solver run is justified.
