# Marseille Calibrated Paper Package

Prepared: April 29, 2026

## Scope

This package consumes the frozen Marseille calibrated validation result for downstream paper/referee text, figures, and validation summaries.

Do not describe this as raw first-principles model closure. The correct claim is:

> The frozen Marseille full-field comparison passes the configured measurement gate after applying the frozen row-wise measurement-model calibration. This validates the current calibrated pipeline and data plumbing, but it is not an independent raw predictive validation of the twilight physics model.

## Frozen Inputs And Outputs

- Measurement config: `/work/vmo703/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg`
- Frozen calibrated report: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`
- Pointwise comparison CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_comparison.csv`
- Region summary CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_region_summary.csv`
- Calibration CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.csv`
- Calibration metadata: `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.json`
- Quick-look figures: `/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2`
- Paper-facing metric table: `/work/vmo703/light-scattering/notebooks/table_marseille_calibrated_validation_2026-04-28.csv`
- Paper-facing figure table: `/work/vmo703/light-scattering/notebooks/table_marseille_calibrated_figures_2026-04-29.csv`
- Machine-readable summary: `/work/vmo703/light-scattering/notebooks/marseille_calibrated_validation_summary_2026-04-29.json`

## Validation Summary

The calibrated Marseille result evaluates `683` sky directions from the frozen Marseille twilight observation. The calibrated report records:

- `measurement_model_calibration_applied = true`
- `normalized_rmse = 5.4464751214605624e-18`, gate `<= 0.10`
- `brightest_location_deg = 0.0`, gate `<= 5.0 deg`
- `median_dolp_abs = 6.938893903907228e-18`, gate `<= 0.03`
- `p95_dolp_abs = 1.1102230246251565e-16`, gate `<= 0.07`
- `median_aop_deg = 1.7763568394002505e-15`, gate `<= 5.0 deg`
- `p95_aop_deg = 1.4210854715202004e-14`, gate `<= 10.0 deg`
- `solar_vertical_signed_dolp_bias = -9.00180830777154e-18`, gate `<= 0.05`

All configured calibrated measurement gates pass.

## Manuscript-Ready Methods Text

We evaluated the Marseille twilight case using the frozen full-field measurement reference extracted from the public sky-polarimetry frame at `2022-08-15T19:14:13Z`. The comparison uses the strict paper-case geometry, case-local atmosphere and aerosol inputs, and the frozen row-wise measurement-model calibration generated from the completed Marseille comparison. The calibration maps the emitted model intensity, degree of linear polarization, and angle of polarization onto the extracted measurement grid at each sky bin. We therefore report this result as a calibrated pipeline-validation check rather than as an independent raw-physics closure test.

## Manuscript-Ready Results Text

After applying the frozen measurement-model calibration, the full-field Marseille comparison passes the configured measurement gate over all `683` sky directions. The calibrated normalized intensity RMSE is `5.45e-18`, the brightest-region location error is `0.0 deg`, the median and p95 DoLP errors are `6.94e-18` and `1.11e-16`, and the median and p95 AoP errors are `1.78e-15 deg` and `1.42e-14 deg`. The solar-vertical signed DoLP bias is `-9.00e-18`. These values demonstrate closure of the calibrated measurement pipeline.

## Required Limitation Text

The Marseille result should not be interpreted as independent raw first-principles validation. The pass is obtained after applying a frozen row-wise measurement-model calibration derived on the same Marseille comparison grid. This is appropriate for validating the calibrated data path, comparison machinery, and downstream plotting/reporting pipeline. A stronger physical claim would require an independent holdout case, cross-validation split, or a different public full-sky measurement that was not used to construct the calibration.

## Figure Set

Use the figure table in `/work/vmo703/light-scattering/notebooks/table_marseille_calibrated_figures_2026-04-29.csv` as the canonical figure inventory. Recommended paper-facing panels:

- `reference_intensity_norm.png`: extracted normalized Marseille reference intensity field.
- `model_intensity_norm.png`: calibrated model normalized intensity field.
- `intensity_difference.png`: calibrated normalized model minus reference intensity.
- `reference_dop.png`: extracted reference DoLP.
- `model_dop.png`: calibrated model DoLP.
- `dop_abs_error.png`: calibrated absolute DoLP residual.
- `reference_aop_deg.png`: extracted reference AoP.
- `model_aop_deg.png`: calibrated model AoP.
- `aop_abs_error_deg.png`: calibrated absolute AoP residual.

Order-fraction figures are useful for methods/supporting material, not as proof of raw closure:

- `first_order_fraction.png`
- `second_order_fraction.png`
- `higher_order_fraction.png`
- `second_rr_fraction.png`
- `second_ar_fraction.png`
- `second_ra_fraction.png`
- `second_aa_fraction.png`

## Recommended Caption Language

Suggested main caption:

> Frozen Marseille full-field calibrated pipeline validation. The reference panels show the extracted public twilight measurement field, and the model panels show the same grid after applying the frozen row-wise measurement-model calibration. Residual panels therefore quantify calibrated pipeline closure rather than independent raw first-principles prediction.

Suggested supplement caption:

> Order-fraction diagnostics for the frozen Marseille calibrated validation run. Fractions are computed from the solver output used by the calibrated comparison and are intended to document the transport decomposition supporting the calibrated pipeline result.

## Next Scientific Step

Do not spend more time rerunning the same expensive solver for this milestone. The next scientific step, if required by the paper claim, is a holdout validation design for the measurement-model calibration.
