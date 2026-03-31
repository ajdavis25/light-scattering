# Monte Carlo Clear-Sky Twilight Solver

This directory now contains the active production solver path for the repository.

Current solver assumptions are documented explicitly in [MODEL_ASSUMPTIONS.md](/c:/Users/ashton/Desktop/projects/light-scattering/MODEL_ASSUMPTIONS.md).

The scientific target is a clear-sky twilight radiance and polarization study with:

- a 1D spherically stratified atmosphere
- observer-centered fisheye sky output
- wavelength-resolved transport
- Stokes `I/Q/U/V` transport through scattering events
- validation hooks for analytic checks, benchmark comparisons, and measurement comparisons

## Current status

The repo is materially beyond the earlier prototype stage. The checked-in validation suite currently passes, but the solver is still not automatically paper-safe for quantitative claims.

What is implemented:

- CSV-driven atmosphere ingestion in [src/Atmosphere.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.hpp) and [src/Atmosphere.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.cpp)
- wavelength grid, solar spectrum ingestion, optional instrument-response weighting, gas absorption tables, and Rayleigh cross sections in [src/WavelengthHandling.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.hpp) and [src/WavelengthHandling.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.cpp)
- Rayleigh and aerosol phase/Mueller handling in [src/PhaseFunctions.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.hpp) and [src/PhaseFunctions.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.cpp)
- a more realistic bundled aerosol reference with a profile-driven scattering/absorption Angstrom split and a generated multi-component clear-sky continental phase table in [data/atmosphere/clear_sky_midlatitude.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/atmosphere/clear_sky_midlatitude.csv), [data/optics/aerosol_phase_matrix_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_phase_matrix_reference.csv), and [data/optics/aerosol_reference_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_reference_metadata.json)
- Stokes basis rotation and Mueller transport in [src/Polarization.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Polarization.hpp) and [src/Polarization.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Polarization.cpp)
- Lambertian land handling plus an optional Cox-Munk ocean surface path in [src/SurfaceReflection.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/SurfaceReflection.hpp) and [src/SurfaceReflection.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/SurfaceReflection.cpp)
- production config, solver entrypoint, and structured sky output in [src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp), [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp), and [src/main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/main.cpp)
- a deterministic single-scatter control variate, survival-biased higher-order transport, and source-guided first higher-order sampling in [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp)
- an optional Rayleigh polarization-guided first-event proposal in [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp), currently left disabled in the default config because the first tested tuning candidate regressed overall convergence
- a deterministic second-order atmospheric control variate in [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp), enabled in the default config with low-order spatial and angular quadrature
- validation runner in [src/ValidationQA.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.hpp), [src/ValidationQA.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.cpp), and [src/validation_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/validation_main.cpp)
- a real frozen twilight paper-case path in [config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg), [frozen_marseille_twilight_20220815_191413z.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg), [frozen_marseille_twilight_20220815_191413z_measurement.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg), and [data/paper_cases/frozen_marseille_twilight_20220815_191413z](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z)
- a paper-case builder in [tools/build_paper_case_inputs.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/build_paper_case_inputs.py) and shared gas cross-section generator in [tools/generate_gas_cross_sections.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_gas_cross_sections.py)
- benchmark/measurement reference loading now supports per-point solar geometry, which is needed for published spherical limb benchmark subsets in [src/ValidationQA.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/ValidationQA.cpp)
- optional OpenMP acceleration is enabled in the current MinGW build through [CMakeLists.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/CMakeLists.txt)
- bundled scalar and vector benchmark cases plus bundled intensity, meridian-polarization, and full-sky-polarization measurement cases now exist under [data/reference_cases](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases)
- two published spherical-vector limb benchmarks from Zawada et al. are now bundled in [config/benchmark_zawada_spherical_vector_single.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_single.cfg) and [config/benchmark_zawada_spherical_vector_multiple.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg), with reference subsets in [data/reference_cases/benchmark_zawada_spherical_vector_single.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_single.csv) and [data/reference_cases/benchmark_zawada_spherical_vector_multiple.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases/benchmark_zawada_spherical_vector_multiple.csv), and shared 351 nm benchmark inputs in [data/atmosphere/benchmark_zawada_limb_351nm.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/atmosphere/benchmark_zawada_limb_351nm.csv), [data/optics/aerosol_phase_matrix_zawada_351nm.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_phase_matrix_zawada_351nm.csv), and [data/optics/ozone_cross_section_zawada_351nm.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/ozone_cross_section_zawada_351nm.csv)

What is still blocking publication use:

- the second-order control-variate path is materially heavier than the previous estimator and can push a full validation run beyond interactive tool timeouts
- the default passing benchmark gate is still the older plane-parallel Rayleigh pair; the stricter published spherical-vector limb benchmarks are now bundled separately, and both Zawada smoke benchmarks now pass their benchmark metrics
- the default passing polarization measurement is still a small meridian twilight subset; a stricter calibrated full-sky fisheye DoLP case is now bundled separately and currently fails
- aerosol optical properties are now more realistic than the earlier placeholders. The default repo tables are still development references, but the frozen Marseille paper case now writes and consumes case-local aerosol optics and a case-local IMX250MYR blue-channel response proxy under [data/paper_cases/frozen_marseille_twilight_20220815_191413z](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z)
- the active modeling assumptions are still constrained to the clear-sky 1D spherical case documented in [../MODEL_ASSUMPTIONS.md](/c:/Users/ashton/Desktop/projects/light-scattering/MODEL_ASSUMPTIONS.md)
- the stricter benchmark story is materially stronger now because the published spherical Zawada single-scatter and all-orders smoke benchmarks both pass their benchmark metrics, but the repo still needs stronger measurement closure before paper use
- the frozen Marseille twilight case now has a real extracted reference field in [data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference.csv), built from the raw public frame and calibration files by [tools/reduce_marseille_case.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/reduce_marseille_case.py)
- that reducer now also consumes the public Marseille `rotation.npy` calibration and rotates each pixel's measured sensor-frame `Q/U` into the comparison basis before binning, using `relative_azimuth_deg - rotation_z_deg - 90 deg` as the per-ray correction. That removes the earlier `45 deg / 215 deg` AoP sector aliasing from the Marseille reference path.
- the Marseille zero-field failure was first traced to an O4 units bug in [tools/generate_gas_cross_sections.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_gas_cross_sections.py) and [data/optics/o4_cross_section_reference.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/o4_cross_section_reference.csv). That bug is fixed, and [tests/test_wavelength.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests/test_wavelength.cpp) now guards against implausible O4 extinction.
- after that O4 fix, the current paper-path blocker is still the Marseille twilight mismatch itself, but the diagnosis is narrower now. This pass added basis-consistent fixed-reference Mueller handling, order-resolved first/second/higher outputs, and a physics-general below-horizon higher-order proposal in [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp).
- the strict paper configs in [config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg), [config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg), and [config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg) now set `strict_paper_mode=true`, and the loader rejects any paper config that still tries to use the old Marseille-only empirical twilight boost fields.
- a Marseille azimuth-convention sweep ruled out a simple `90 deg` frame mismatch as the dominant issue. The unshifted reference remains the best match; rotated variants are worse.
- step 7 is now implemented in the Marseille paper-case builder. It fetches matched Open-Meteo pressure-level thermodynamics, uses Marseille_ATMO AERONET direct-sun AOD for the column target, and upgrades the case-local aerosol model to an inversion-backed fallback when no strict within-window AERONET inversion exists. The resulting [data/paper_cases/frozen_marseille_twilight_20220815_191413z/aerosol_optics.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/aerosol_optics.json) now records `model_name = marseille_case_local_aeronet_inversion_v2`, `aeronet_inversion_usage_mode = same_day_fallback`, and the inversion-derived `SSA`, asymmetry, refractive-index, depolarization, and fine/coarse fractions at `550 nm`
- the Marseille instrument response is now tracked as a public-doc-constrained IMX250MYR proxy rather than an unqualified internal curve. [data/paper_cases/frozen_marseille_twilight_20220815_191413z/instrument_response.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/instrument_response.csv), [data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_metadata.json), and [data/paper_cases/frozen_marseille_twilight_20220815_191413z/paper_case_provenance.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/paper_case_provenance.json) now cite the Sony IMX250 Polarsens flyer and the LUCID PHX050S polarized technical reference as the public basis for that proxy
- the Marseille reference reduction now writes `q`, `u`, `DoLP`, and `AoP`, and the measurement comparison outputs now include AoP errors, signed DoLP bias, and first/second/higher order fractions
- the solver now exposes explicit runtime-budget knobs for deterministic integrals, and the new full-field strict interactive config [config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_interactive.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_interactive.cfg) completes on the full `683`-bin Marseille field for interactive debugging
- the remaining Marseille blocker is still mainly polarization closure in below-horizon twilight transport; the strict paper path no longer depends on the empirical Marseille boost, but a full strict Marseille batch rerun still needs to be completed outside the interactive shell budget

## Directory layout

- [config](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config): reproducible simulation configs
- [data/atmosphere](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/atmosphere): atmosphere profile CSVs
- [data/optics](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics): solar spectrum, ozone, aerosol, and surface tables
- [data/reference_cases](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/reference_cases): external benchmark and measurement inputs
- [data/paper_cases](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases): frozen or interim paper-case packages
- [tools/generate_reference_cases.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_reference_cases.py): reproducible generator for the bundled benchmark and measurement CSV/metadata files
- [tools/generate_aerosol_reference.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_aerosol_reference.py): reproducible generator for the bundled aerosol profile and aerosol phase/Mueller lookup table
- [tools/build_paper_case_inputs.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/build_paper_case_inputs.py): builds the frozen Marseille twilight paper case from matched public inputs and writes the provenance bundle
- [tools/reduce_marseille_case.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/reduce_marseille_case.py): downloads the raw Marseille day file plus calibration assets, extracts the exact frozen frame, and reduces it into the validator CSV format
- [tools/generate_gas_cross_sections.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_gas_cross_sections.py): writes shared O2/O4/H2O/NO2 reference cross-section tables
- [tools/run_paper_validation_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_paper_validation_batch.py): launches the strict paper-validation stack with a persistent log, including a detached background mode for the full Marseille-heavy batch run
- [tools/run_marseille_paper_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_marseille_paper_batch.py): runs the strict frozen Marseille measurement case first and then [config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg) under one persistent batch log
- [results](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results): generated solver and validation outputs
- [tests](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tests): executable physics/unit-style checks

## Build

From the repository root on Windows with MinGW:

```powershell
$env:Path = "C:\mingw64\bin;" + $env:Path
& "C:\Program Files\CMake\bin\cmake.exe" -S monte_carlo_cpp -B monte_carlo_cpp\build_current -G "MinGW Makefiles" `
  -DCMAKE_MAKE_PROGRAM=C:/mingw64/bin/mingw32-make.exe `
  -DCMAKE_C_COMPILER=C:/mingw64/bin/gcc.exe `
  -DCMAKE_CXX_COMPILER=C:/mingw64/bin/g++.exe
& "C:\Program Files\CMake\bin\cmake.exe" --build monte_carlo_cpp\build_current
```

If the local MinGW driver path is failing during the assembler or linker handoff, use the checked-in workaround builder instead:

```powershell
python monte_carlo_cpp\tools\build_with_mingw_workaround.py
```

That script rebuilds [build_current](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/build_current) by compiling to assembly with `g++`, assembling with `as.exe`, linking with `ld.exe`, and staging the required MinGW runtime DLLs next to the executables.

To regenerate the frozen paper case and shared gas tables:

```powershell
python monte_carlo_cpp\tools\generate_gas_cross_sections.py
python monte_carlo_cpp\tools\build_paper_case_inputs.py
```

To launch the strict paper-validation stack as a logged batch job:

```powershell
python monte_carlo_cpp\tools\run_paper_validation_batch.py --detached
```

That writes a persistent log to [results/validation/paper_validation_batch.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/paper_validation_batch.log) and leaves the final machine-readable metrics in [results/validation/validation_report.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/validation_report.json).

## Run the production solver

```powershell
.\monte_carlo_cpp\build_current\MonteCarloCPP.exe monte_carlo_cpp\config\default_clear_sky.cfg
```

That writes:

- [results/default_clear_sky/sky_result.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky/sky_result.csv)
- [results/default_clear_sky/sky_result_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/default_clear_sky/sky_result_metadata.json)

## Run validation

```powershell
.\monte_carlo_cpp\build_current\ValidationRunner.exe monte_carlo_cpp\config\default_clear_sky.cfg
```

That writes:

- [results/validation/validation_report.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/validation_report.json)

Current validation result:

- analytic Rayleigh checks pass
- unit-style tests pass
- convergence currently passes with `convergence_peak_intensity_rel = 0.0382187`, `convergence_peak_dolp_abs = 0.00138362`, and `convergence_flux_rel = 0.00211343`
- the scalar and vector benchmark comparisons now both run and pass
- the intensity and polarization measurement comparisons now both run and pass
- the stricter published spherical-vector limb benchmarks are now bundled in [config/benchmark_zawada_spherical_vector_single.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_single.cfg) and [config/benchmark_zawada_spherical_vector_multiple.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg), and both are included in [config/default_clear_sky_extended_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/default_clear_sky_extended_validation.cfg)
- the stronger all-orders spherical case now also has a reproducible smoke wrapper in [config/benchmark_zawada_spherical_vector_multiple_smoke.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple_smoke.cfg)
- the benchmark-specific aerosol extinction generation bug is now fixed in [tools/generate_reference_cases.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/generate_reference_cases.py): `cm^-1` is now converted to `m^-1` with `1e2`, not `1e-2`
- the bounded smoke run for the single-scatter spherical case now passes with `median_intensity_rel = 0.000459695`, `p95_intensity_rel = 0.00140769`, `median_dolp_abs = 3.37895e-05`, and `p95_dolp_abs = 0.000292954` in [benchmark_zawada_spherical_vector_single_smoke.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/benchmark_zawada_spherical_vector_single_smoke.txt)
- the bounded smoke run for the stronger all-orders spherical case now also passes its benchmark metrics with `median_intensity_rel = 0.00457567`, `p95_intensity_rel = 0.012131`, `median_dolp_abs = 0.0199496`, and `p95_dolp_abs = 0.0353911` in [benchmark_zawada_spherical_vector_multiple_smoke.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/benchmark_zawada_spherical_vector_multiple_smoke.txt)
- the stronger all-orders benchmark now uses a benchmark-only heavier configuration in [benchmark_zawada_spherical_vector_multiple.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple.cfg): deterministic first/second-order pullout, denser second-order quadrature, and modest Rayleigh polarization-guided first-event sampling
- the stricter full-sky fisheye polarization case is not part of the default passing gate; it is evaluated separately and currently fails with `median_dolp_abs = 0.0660017` and `p95_dolp_abs = 0.217099`
- the exact strict Gal-case rerun now works again from the rebuilt workaround-based [build_current](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/build_current), and the pointwise diagnostics are now written to [measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv)
- the paper-path validation config now exists in [config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg), but it is expected to fail until the primary twilight full-sky dataset is frozen and the full batch run completes
- the full default validation stack remains heavy enough to exceed the interactive shell timeout in the current environment, so the strict end-to-end rerun is still best treated as batch work rather than an interactive command

## Measurement-only case evaluation

For stricter published measurement checks without paying for the entire validation stack:

```powershell
.\monte_carlo_cpp\build_current\MeasurementCaseRunner.exe monte_carlo_cpp\config\measurement_gal_lapland_fullsky_450nm_dolp.cfg
```

Current exact measurement-case reports are stored in:

- [results/measurement_case_reports/measurement_koomen_meridian_hminus6_polarization.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_koomen_meridian_hminus6_polarization.txt)
- [results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp.txt)
- [results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv)
- [results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive.txt)
- [results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_comparison.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_comparison.csv)
- [results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv)

The bundled strict full-sky fisheye case is configured in [config/measurement_gal_lapland_fullsky_450nm_dolp.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/measurement_gal_lapland_fullsky_450nm_dolp.cfg). An extended default validation config that includes it is available at [config/default_clear_sky_extended_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/default_clear_sky_extended_validation.cfg), but that stricter path is expected to fail with the current solver/input stack.

The current exact Gal-case report is:

- `median_dolp_abs = 0.0660017`
- `p95_dolp_abs = 0.217099`
- `max_dolp_abs = 0.308684` at `zenith_deg = 45`, `relative_azimuth_deg = 270`

Quick-look measurement comparison plots for that case are now written by the Python layer under [../plots/current/measurement_cases/measurement_gal_lapland_fullsky_450nm_dolp](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/measurement_cases/measurement_gal_lapland_fullsky_450nm_dolp).

For Marseille-specific solver work, the reducer now also writes a strict-physics subset reference at [data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_strict_subset.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_strict_subset.csv) and a matching config at [config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_strict_subset.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_strict_subset.cfg). That subset keeps the strict paper physics and full spectral band, but reduces the scored directions and Monte Carlo budget for iteration.

The reducer now also writes a strict-paper profiling subset at [data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_profile_subset.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_profile_subset.csv) with summary metadata in [data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_profile_subset_summary.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference_profile_subset_summary.json) and a matching config at [config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset.cfg). This is the right runtime-profiling surface because it keeps the strict paper physics and only reduces the number of scored Marseille directions.

The measurement runner now also writes region-level Marseille diagnostics. On the current interactive Marseille run, the strongest remaining miss is still the solar-vertical mid-zenith band, and that is visible directly in [frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_region_summary.csv).

[src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp) now emits live runtime instrumentation for measurement runs: a startup line, periodic per-direction progress lines, and timing summary fields in the final report. A profiling run of the strict-paper Marseille subset writes to [results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset_progress.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset_progress.log). That log now emits completed first-order and second-order stage lines on the strict subset, so the Marseille paper path is still expensive but no longer opaque.

The profiler now instruments the direction solve internally as well. [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) emits live first-order band progress, second-order band progress, and second-order incoming-ray timing through the measurement callback. The current Marseille strict-profile log shows that the first runtime wall is deterministic single scatter before second order becomes visible: the first band takes about `3.9-5.3 s` for near-zenith directions, about `11.7-13.3 s` for `~69-83 deg` zenith directions, and about `30.6 s` for the `87.65 deg` horizon-skimming direction. The same log exposes the underlying driver: `first_order_steps` ranges from `101` up to `902`, so the strict Marseille cost currently grows first with line-integral depth and solar-disk sampling before deterministic second scatter becomes the leading bottleneck.

The first deterministic single-scatter runtime reduction is now checked in as well. [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) now caches solar-disk nodes in the simulation context, prunes spectrally negligible instrument-weighted bands from the active Marseille paper-case band set, precomputes sun-ray optical-depth spectra once per view-step and solar-disk node rather than rebuilding the same sun-transmittance geometry separately for every band, and then integrates those sun-path optical depths by spherical altitude-shell crossings instead of a fixed-distance brute-force path walk. On the current Marseille blue-channel case, the strict active spectral set stays at `25` bands, but the measured runtime shape changes by roughly an order of magnitude: in [results/measurement_case_reports/profile_subset_runtime_probe.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/profile_subset_runtime_probe.log), the near-zenith `101-step` directions now complete first order in about `0.90-0.94 s`, the `117-step` directions in about `0.64-0.97 s`, the `263-step` directions in about `1.12-2.12 s`, a `595-step` direction in about `2.31 s`, and the `902-step` horizon-skimming direction in about `4.82 s`. That means the low-elevation deterministic first-order wall is no longer the dominant strict Marseille bottleneck; the remaining runtime pressure has shifted to deterministic second order on the deepest low-elevation directions.

The next deterministic runtime reduction is now checked in too. [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) no longer solves the incoming deterministic single-scatter ray separately for every spectral band inside the second-order quadrature. It now computes the full active-band incoming single-scatter spectrum once per `(view step, mu, phi)` geometry and reuses that cached result across the whole band loop. On the same strict Marseille profiling subset, that cuts the second-order wall heavily: the current [results/measurement_case_reports/profile_subset_runtime_probe.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/profile_subset_runtime_probe.log) now shows `second_order_incoming_single_scatter_calls = 72` for the `6x4x3x4` cases, `196` for the `7x5x4x7` cases, and `320` for the `8x6x5x8` cases, with `second_inner_s` dropping to about `0.62-0.78 s` for the low-zenith directions, about `3.13-3.57 s` for the `~69 deg` directions, and about `5.20-5.84 s` for the deepest `~83-88 deg` directions. The strict Marseille subset is therefore no longer suffering from the original first/second-order runtime explosion.

The earlier Marseille AoP quadrant/sign failure in the `45 deg` and `215 deg` relative-azimuth sectors is materially reduced. After the calibrated basis rotation was added to the reducer, the interactive full-field report now gives `region_aop_flip_relaz_45_sector_median_aop_deg = 16.71` and `region_aop_flip_relaz_215_sector_median_aop_deg = 16.3629`, rather than the earlier `~80-90 deg` aliasing pattern.

The deterministic second-order path now also has a physics-general below-horizon adaptive quadrature. In twilight, [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) can raise the effective second-order `view_steps`, `ray_steps`, `mu_nodes`, and `phi_nodes` above the config base values without relying on any Marseille-only empirical boost. The Marseille measurement runner exposes the result through additional region metrics: `mean_second_rr_frac`, `mean_second_ar_frac`, `mean_second_ra_frac`, and `mean_second_aa_frac`.

With that step-4 upgrade, the interactive full-field Marseille run is no longer first-order-only in the worst regions. The current report now shows `region_solar_vertical_midzen_mean_first_frac = 0.535033`, `region_solar_vertical_midzen_mean_second_frac = 0.464967`, `region_bright_horizon_arc_mean_first_frac = 0.472768`, and `region_bright_horizon_arc_mean_second_frac = 0.527232`. The overall interactive Marseille metrics also improved to `normalized_rmse = 0.119553` and `median_dolp_abs = 0.163825`, though that is still far from the paper gate.

Step 5 is now implemented in the solver as well. The third-and-higher continuation guide is no longer accidentally tied to the first-scatter source-guided gate, and [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) now uses stratified phase/tangent/horizon continuation branches with MIS-consistent mixture pdfs for the third-plus twilight remainder. The Marseille paper configs use `twilight_higher_order_branches = 4`, while the interactive full-field config keeps `twilight_higher_order_branches = 3` so it remains runnable locally.

On the current interactive Marseille rerun, that step-5 path makes the higher-order fraction measurably nonzero in the key biased regions: `region_solar_vertical_midzen_mean_higher_frac = 0.0126252`, `region_bright_horizon_arc_mean_higher_frac = 0.0550218`, and `region_antisolar_midzen_mean_higher_frac = 0.0634904`. The overall interactive metrics after the step-5 pass are `normalized_rmse = 0.155497`, `median_dolp_abs = 0.169169`, and `solar_vertical_signed_dolp_bias = 0.374188`. That is still not paper-safe, but it confirms the remaining Marseille blocker is not a missing third-plus continuation path anymore.

Step 6 is now implemented as a reproducible surface sensitivity check. [tools/run_surface_sensitivity.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_surface_sensitivity.py) runs the Marseille measurement config twice, once with `lambertian_land` and once with `coxmunk_ocean`, and writes the result to [results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.txt) and [results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_interactive_surface_sensitivity_summary.json).

The current Marseille result is `max_abs_delta = 0.0`, so surface choice is negligible on the present twilight path. That is not evidence that surface physics can never matter; it reflects the current model structure. In [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp), surface radiance enters only through the direct-solar ground reflection path. For the frozen Marseille case, the Sun is below the local horizon, so switching between land and Cox-Munk ocean does not currently change the modeled field.

After the step-7 input upgrade, the frozen Marseille case artifacts on disk have been regenerated with the inversion-backed aerosol model and updated provenance. A fresh full-field Marseille rerun under those inputs still needs a controlled batch execution, because interactive `MeasurementCaseRunner.exe` sessions continue to hit shell time limits on this machine before exiting cleanly.

The repo now includes a resumable exact-direction batch runner at [tools/run_measurement_case_batched.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py). It splits a measurement reference into exact-direction chunks, preserves the original direction indices for reproducible Monte Carlo seeding, and writes checkpointed partial outputs under `results/measurement_case_reports`. The Marseille paper-batch launcher in [tools/run_marseille_paper_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_marseille_paper_batch.py) now uses that batch path for the measurement phase.

That exact-direction path now also supports per-direction checkpointing inside a single Marseille direction. [src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp) can now reuse completed deterministic first/second-order results and resume only the third-plus Monte Carlo remainder, while [src/measurement_case_main.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp) exposes that through `--checkpoint-state` and `--higher-order-block-size`. The current batch runner uses this automatically when `--batch-size 1`.

That checkpoint/resume path is already verified on the lightweight Marseille diagnostic config [config/paper_cases/fast16_base.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/fast16_base.cfg). With `batch_size = 1` and `higher_order_block_size = 1`, the first invocation computed deterministic first/second order plus one higher-order sample in about `38.38 s`, and the second invocation resumed from the saved checkpoint and finished the remaining higher-order sample in about `3.37 s`. The merged row is recorded in [results/measurement_case_reports/fast16_base_batched_partial_rows.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/fast16_base_batched_partial_rows.csv).

The batched log path is now trustworthy as a live progress surface as well. [tools/run_measurement_case_batched.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py) now line-flushes child output into the batched log, and [tools/run_marseille_paper_batch.py](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/tools/run_marseille_paper_batch.py) now launches the measurement runner under unbuffered Python. The current [results/measurement_case_reports/fast16_base_batched.log](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/fast16_base_batched.log) shows the checkpoint progress lines directly, so `Get-Content -Wait` can now be used as a real progress monitor for the batched measurement path.

## Run tests

```powershell
& "C:\Program Files\CMake\bin\ctest.exe" --test-dir monte_carlo_cpp\build_current --output-on-failure
```

## Python plotting/comparison layer

The Python layer is now downstream only. It reads production outputs, writes an NPZ bundle, and generates comparison plots against the analytic single-scatter reference.

Run it from the repository root:

```powershell
python -m spherical.main
```

Current plots are written under [..\plots\current](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current).

The default measurement plotting target in [..\spherical\main.py](/c:/Users/ashton/Desktop/projects/light-scattering/spherical/main.py) is now the full-field interactive Marseille case, not the older Gal case. The Marseille measurement plot set now includes DoLP, AoP, signed DoLP bias, normalized intensity, and order-fraction fields under [..\plots\current\measurement_cases](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/measurement_cases).

## Production config interface

The active config schema is centered on the `SimulationConfig` interface in [src/MonteCarloDriver.hpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.hpp). The default checked-in config is [config/default_clear_sky.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/default_clear_sky.cfg).

The important configuration groups are:

- observer location and altitude
- solar geometry
- atmosphere and optics input tables
- spectral grid
- surface albedo table
- Monte Carlo controls such as photon count, random seed, and Russian roulette threshold
- a `deterministic_single_scatter` switch that keeps the single-scatter baseline deterministic and leaves Monte Carlo to estimate the higher-order remainder
- source-guided first-scatter controls for the higher-order twilight estimator
- output grid resolution
- benchmark and measurement case config paths for validation
- validation thresholds

## Paper gate

Do not use outputs from this directory as publication figures or quantitative paper tables until all of the following are true:

1. the convergence metrics in [results/validation/validation_report.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/validation/validation_report.json) pass
2. the bundled or study-specific external benchmark case passes the benchmark gate
   For paper claims, this should include the stricter spherical-vector benchmark in [config/default_clear_sky_extended_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/default_clear_sky_extended_validation.cfg), not only the default plane-parallel benchmark pair.
3. at least one documented public/reference twilight measurement case passes the measurement gate
4. the current bundled benchmark and measurement cases are replaced or supplemented with study-appropriate validation data for the specific paper claim
5. the atmosphere and aerosol optical inputs are replaced with defensible study-specific datasets
6. [config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg) passes with `paper_primary_measurement_frozen = true` and the primary twilight full-sky measurement case promoted out of interim status
