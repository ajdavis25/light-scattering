# Model Assumptions

Last updated: 2026-04-29

This file records the assumptions that are currently hard-coded or implied by the active production path in:

- [monte_carlo_cpp/src/MonteCarloDriver.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp)
- [monte_carlo_cpp/src/Atmosphere.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Atmosphere.cpp)
- [monte_carlo_cpp/src/WavelengthHandling.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/WavelengthHandling.cpp)
- [monte_carlo_cpp/src/PhaseFunctions.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/PhaseFunctions.cpp)
- [monte_carlo_cpp/src/Polarization.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/Polarization.cpp)
- [monte_carlo_cpp/src/SurfaceReflection.cpp](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/src/SurfaceReflection.cpp)

The assumptions below describe the current solver, not the final research target.

## April 29, 2026 Marseille Calibration Addendum

The frozen Marseille full-field result now passes the configured measurement gate only after applying a row-wise measurement-model calibration:

- report: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`
- calibration: `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.csv`
- paper package: `/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md`

This is a calibrated pipeline-validation result. It should not be cited as independent raw first-principles closure of the below-horizon twilight polarization model.

## Physical Assumptions

### Atmosphere

- The atmosphere is **clear sky only**.
- The atmosphere is **1D, spherically stratified, and horizontally homogeneous**.
- Atmospheric state is read from a single vertical profile CSV.
- Positive profile quantities such as pressure, molecular density, ozone density, and aerosol extinction are interpolated in log space with altitude.
- Temperature, single-scattering albedo, asymmetry, and aerosol spectral exponents are interpolated linearly with altitude.
- The Earth is a sphere with fixed radius `6371 km`.
- The top of atmosphere is taken from the profile top or the configured override, with `100 km` used in the default clear-sky case.
- There is **no 3D structure**, no local aerosol plumes, no clouds, no terrain-shadow field, and no horizontal gradients.

### Optical Constituents

- Molecular scattering is Rayleigh only.
- The default development path still uses ozone as the only enabled gas-absorption table.
- The paper-path configs can now also load O2, O4, H2O, and NO2 absorption tables.
- Aerosol extinction is taken from the atmosphere profile at `550 nm`.
- Aerosol spectral scaling is split into profile-driven scattering and absorption Angstrom exponents.
- Aerosol single-scattering albedo and asymmetry are profile-driven.
- The default bundled aerosol profile now represents a multi-component clear-sky continental mixture with boundary-layer fine and coarse modes, an elevated fine mode, and a weak background sulfate layer.
- Aerosol scattering matrix elements come from the checked-in aerosol phase/Mueller lookup table documented in [monte_carlo_cpp/data/optics/aerosol_reference_metadata.json](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/optics/aerosol_reference_metadata.json).

### Solar Source

- The source is the **direct Sun only**.
- The default development path treats the Sun as a directional beam.
- The paper-path configs can now enable a finite solar disk using a fixed 7-node equal-area quadrature with no limb darkening.
- There is no moonlight, starlight, airglow, artificial light, or surface thermal emission.
- The repo now also contains an explicit strict interactive Marseille config that keeps the full field but lowers the deterministic integral and spectral budgets for debugging. That config is for engineering turnaround only, not paper metrics.
- Solar geometry is either:
  - explicit solar zenith and azimuth, or
  - observer latitude plus solar declination and hour angle.
- No external ephemeris package is currently used in the production C++ path.

### Surface

- The ground is a **Lambertian reflector**.
- Surface reflectance is spectral but spatially uniform within a run.
- The default production run uses the checked-in land albedo table.
- The paper-path configs now support `lambertian_land` and `coxmunk_ocean` surface models.
- In the current twilight transport implementation, surface radiance enters through the direct-solar ground-reflection path only.
- For below-horizon twilight cases like the frozen Marseille case, that means the current land-versus-ocean choice is effectively inactive unless the surface transport model is extended beyond direct solar reflection.
- There is still no terrain slope model or subpixel land-cover mixture.

### Polarization

- The active solver transports Stokes `I/Q/U/V`.
- Linear polarization is modeled through Mueller matrices and basis rotations at scattering events.
- Circular polarization `V` is tracked numerically but is not expected to be important in the intended clear-sky cases.
- The passing default validation suite includes both a vector polarization benchmark and a twilight meridian polarization measurement subset.
- A stricter full-sky calibrated polarization measurement case is now bundled, but the current solver does not yet match it closely enough for it to be part of the default passing gate.

### Geometry and Output

- Output is observer-centered over the **upper hemisphere** only.
- Default output grid is a uniform zenith/azimuth binning of `19 x 36`.
- Zenith angle is measured from the local zenith:
  - `0 deg` = overhead
  - `90 deg` = horizon
- Azimuth is measured in the local horizontal frame with:
  - `0 deg` = north
  - `90 deg` = east
  - `180 deg` = south
  - `270 deg` = west
- The fisheye rendering pipeline is downstream of the C++ result files; the C++ solver is the scientific source of truth.

## Numerical / Estimator Assumptions

### Monte Carlo Formulation

- The active solver is **backward Monte Carlo** from observer sky bins.
- The direct single-scatter contribution is currently handled with a deterministic control variate.
- Higher-order transport is currently estimated stochastically.
- The first higher-order continuation step now uses:
  - source-guided proposal sampling
  - a phase/sun-cone mixture proposal
  - branch splitting at the first higher-order event
- Russian roulette is used at low throughput with the current fixed threshold from config.
- Survival biasing is used for higher-order scattering instead of analog absorption termination at every real collision.

### Optical Depth and Path Integration

- Optical depth to the Sun is evaluated by numerical line integration using fixed spatial substeps.
- Deterministic single-scatter accumulation along the observer path is also evaluated by numerical substepping.
- These path integrals are not currently adaptive in step size.
- Atmospheric properties between profile levels are interpolated with the same mixed log/linear rules used by the production atmosphere loader.

### Spectral Treatment

- The default production run uses `350-800 nm` with `10 nm` spacing.
- Solar irradiance is read from the checked-in reference spectrum CSV.
- The paper-path configs can now apply an explicit instrument-response CSV when building band weights.
- The final sky result is the sum over the configured spectral bands.
- The bundled output is still not validated against a frozen instrument response from a final paper-grade twilight dataset.

### Validation Gates

- Convergence is currently tested by comparing:
  - a coarse run with half the default photons per bin
  - a fine run with the default photons per bin
- Current thresholds in the validation suite are:
  - `convergence_peak_intensity_rel <= 0.05`
  - `convergence_peak_dolp_abs <= 0.02`
  - `convergence_flux_rel <= 0.02`
- Benchmark and measurement validation currently use normalized field comparisons, not a full absolute radiometric closure study.
- The frozen Marseille full-field pass uses a row-wise measurement-model calibration. That gate validates calibrated pipeline closure, not raw predictive closure.

## Validation Coverage Assumptions

### External Benchmark

- The default passing benchmark suite includes:
  - a scalar principal-plane Rayleigh case generated with PythonicDISORT
  - a vector Rayleigh case from the IPRT A1 intercomparison
- Those benchmark cases validate transport normalization and vector polarization handling under controlled assumptions.
- The bundled vector benchmark is still a **plane-parallel Rayleigh** benchmark, not a full spherical twilight multiple-scattering benchmark.

### Public Measurement Cases

- The default passing measurement suite includes:
  - the Rozenberg (1952) twilight intensity pattern
  - the Koomen et al. (1952) meridian twilight polarization subset
- An additional stricter measurement case is bundled:
  - a coarse full-sky 450 nm DoLP field digitized from Gal et al. (2001) Figure 2a using calibrated full-sky imaging polarimetry
- That Gal et al. case is a full-sky fisheye polarization reference, but it is low-sun daytime sky at `solar_zenith_deg = 83.1`, not below-horizon twilight.
- The current solver does not yet pass that full-sky fisheye polarization case.
- The frozen Marseille full-sky twilight reference is available as a public-data paper case and now passes only as calibrated pipeline validation after applying the frozen row-wise measurement-model calibration.

## What Is Still Missing Before Paper-Safe Use

The current solver is not paper-safe yet because of the following gaps:

- The default passing external benchmark story still tops out at plane-parallel Rayleigh validation for vector polarization, even though stricter published spherical Zawada smoke benchmarks are now bundled separately.
- The stricter bundled full-sky fisheye polarization case currently fails with `median_dolp_abs = 0.0660017` and `p95_dolp_abs = 0.217099` in [monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp.txt](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp.txt).
- Exact pointwise diagnostics for that case now exist in [measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/results/measurement_case_reports/measurement_gal_lapland_fullsky_450nm_dolp_comparison.csv), and quick-look field plots now exist under [plots/current/measurement_cases/measurement_gal_lapland_fullsky_450nm_dolp](/c:/Users/ashton/Desktop/projects/light-scattering/plots/current/measurement_cases/measurement_gal_lapland_fullsky_450nm_dolp).
- The frozen Marseille calibrated plot set now exists under [plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2](/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2).
- The default aerosol optical tables are more realistic than the earlier placeholders, but they are still checked-in reference inputs rather than study-specific externally defended inputs.
- Absolute radiometric credibility still depends on validating the current input tables and spectral handling against the intended paper use case.
- The strongest bundled calibrated full-sky benchmark-style case is still digitized from a published map, but a real raw-public-instrument Marseille twilight case now exists separately under the paper-case workflow.
- A real frozen paper-case package now exists under [monte_carlo_cpp/config/paper_validation.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_validation.cfg), [frozen_marseille_twilight_20220815_191413z.cfg](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z.cfg), and [monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z](/c:/Users/ashton/Desktop/projects/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z).
- That frozen case now uses case-local Marseille aerosol optics built from matched public inputs, including Marseille_ATMO AERONET direct-sun AOD and a same-day AERONET inversion fallback when no strict within-window inversion exists. The case also uses a matched Open-Meteo pressure-level thermodynamic profile with a template upper tail and a public-doc-constrained IMX250MYR blue-channel response proxy. The current paper-facing pass is calibrated; independent raw polarization closure still requires holdout validation or a raw-model improvement.
- The strict paper configs now explicitly reject the Marseille-specific empirical twilight boost path. Those empirical boost knobs remain solver diagnostics only and are not valid in `strict_paper_mode=true`.

## Best Current Interpretation

If the paper is framed as a **development-stage clear-sky twilight radiative transfer solver with calibrated Marseille pipeline validation**, the repo now has a defensible methods/result trail.

If the paper is framed as a **quantitative atmospheric optics study using the simulated twilight fields as independent raw scientific data**, the repo still needs:

1. A stronger vector/multiple-scattering benchmark story than the current plane-parallel Rayleigh benchmark coverage.
2. Holdout or independent validation of the Marseille row-wise calibration.
3. Study-specific aerosol and atmosphere inputs.
4. Absolute/unit-aware validation appropriate to the actual paper claim.
