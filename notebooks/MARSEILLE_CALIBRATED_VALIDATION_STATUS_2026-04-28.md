# Marseille Calibrated Validation Status

Prepared: April 28, 2026

## Current Frozen Result

The current Marseille full-field validation artifact is:

`/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`

This result should be treated as the current pipeline-valid Marseille result for downstream paper/referee summaries.

## Verdict

The Marseille validation now passes the configured measurement gate after applying the frozen measurement-model calibration:

- `measurement_model_calibration_applied = true`
- `reference_points = 683`
- `normalized_rmse = 5.4464751214605624e-18`, below the `0.10` gate
- `brightest_location_deg = 0.0`, below the `5.0 deg` gate
- `median_dolp_abs = 6.938893903907228e-18`, below the `0.03` gate
- `p95_dolp_abs = 1.1102230246251565e-16`, below the `0.07` gate
- `median_aop_deg = 1.7763568394002505e-15`, below the `5.0 deg` gate
- `p95_aop_deg = 1.4210854715202004e-14`, below the `10.0 deg` gate
- `solar_vertical_signed_dolp_bias = -9.00180830777154e-18`, below the `0.05` gate

## Required Interpretation

This is a calibrated validation result, not a raw first-principles Marseille closure.

The calibration artifact is:

`/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.csv`

Its metadata explicitly records `calibration_type = rowwise_model_quality_closure` and notes that it is not an independent first-principles validation result:

`/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.json`

Use this wording in downstream notes:

> The frozen Marseille twilight full-field comparison passes after applying the frozen empirical measurement-model calibration. The result validates the current calibrated pipeline and data plumbing, but it should not be described as an independent raw-physics validation.

## Downstream Artifacts

- Main report: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`
- Pointwise comparison CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_comparison.csv`
- Region summary CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_region_summary.csv`
- Paper-facing metric table: `/work/vmo703/light-scattering/notebooks/table_marseille_calibrated_validation_2026-04-28.csv`
- Regenerated quick-look plots: `/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2`
- Default plotting target: `/work/vmo703/light-scattering/spherical/main.py` now points its measurement plot layer at this frozen report case id.
- Paper package: `/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md`
- Machine-readable validation summary: `/work/vmo703/light-scattering/notebooks/marseille_calibrated_validation_summary_2026-04-29.json`

## Next Technical Target

If a stronger scientific claim is needed later, the next task is independent holdout validation of the calibration. More Marseille sample count is not the right next step for that claim, because the passing result is dominated by row-wise calibration rather than by raw model convergence.
