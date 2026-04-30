# Marseille Light-Scattering Handoff

Prepared: April 30, 2026

## Current Milestone

The current Marseille milestone is packaged as calibrated pipeline validation. Use the frozen full-field Marseille result as the paper-facing artifact only at that claim level:

> The Marseille full-field comparison passes as calibrated pipeline validation after applying the frozen row-wise measurement-model calibration; it should not be described as independent raw first-principles closure.

The result validates the calibrated data path, comparison-grid plumbing, metric extraction, plotting, and configured gate passage. It does not establish independent raw predictive closure of the twilight physics model.

## Canonical Artifacts

- Status note: `/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_VALIDATION_STATUS_2026-04-28.md`
- Paper package: `/work/vmo703/light-scattering/notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md`
- Model assumptions: `/work/vmo703/light-scattering/MODEL_ASSUMPTIONS.md`
- Frozen report: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`
- Pointwise comparison CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_comparison.csv`
- Region summary CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_region_summary.csv`
- Calibration CSV: `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.csv`
- Calibration metadata: `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.json`
- Quick-look figure directory: `/work/vmo703/light-scattering/plots/current/measurement_cases/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2`
- Main paper/referee figure: `/work/vmo703/light-scattering/notebooks/figures/marseille_calibrated_main_panel_2026-04-30.png`
- Machine-readable summary: `/work/vmo703/light-scattering/notebooks/marseille_calibrated_validation_summary_2026-04-29.json`

## Frozen Metrics

Use the full-precision values in the frozen report or machine-readable summary for tables. The current headline values are:

| Quantity | Value |
| --- | ---: |
| Reference sky directions | 683 |
| Measurement-model calibration applied | true |
| Normalized intensity RMSE | 5.4464751214605624e-18 |
| Brightest-location error | 0.0 deg |
| Median DoLP absolute error | 6.938893903907228e-18 |
| p95 DoLP absolute error | 1.1102230246251565e-16 |
| Median AoP error | 1.7763568394002505e-15 deg |
| p95 AoP error | 1.4210854715202004e-14 deg |

All configured calibrated measurement gates pass after applying the frozen row-wise measurement-model calibration.

## Immediate Work Policy

1. Treat the calibrated Marseille package as complete for the current milestone.
2. Use the frozen report, comparison CSV, region summary CSV, calibration CSV, and calibration metadata as the reproducibility bundle.
3. Use the paper package for manuscript-ready wording, referee-response language, limitation text, and figure inventory.
4. Use the main composite figure and the quick-look figure directory listed above for paper/referee packaging.
5. Do not rerun the expensive full Marseille solver unless the calibration, raw model, input package, or paper claim changes.

## Future Scientific Work

If the paper claim remains calibrated pipeline validation, continue with writing and packaging from the frozen artifacts.

If a stronger raw or predictive claim is required, the next scientific task is holdout validation rather than another identical Marseille rerun:

1. Design a withheld-direction calibration split or identify an independent public full-sky twilight measurement.
2. Build and evaluate the calibrated model on sky directions or measurements not used to construct the row-wise calibration.
3. Decide whether additional raw-solver development is justified only after the holdout result is known.
4. Keep historical raw Marseille mismatches as diagnostic and limitation material, especially for below-horizon polarization closure.

## Verification When Code Or Data Changes

For future code or data changes, run the lightest relevant checks first:

```bash
ctest --test-dir monte_carlo_cpp/build_current --output-on-failure
```

For default validation when runtime allows:

```bash
./monte_carlo_cpp/build_current/ValidationRunner monte_carlo_cpp/config/default_clear_sky.cfg
```

For downstream plots:

```bash
python -m spherical.main
```

Avoid full strict Marseille reruns in an interactive shell. If a strict rerun is genuinely needed, use the batch tooling:

```bash
python monte_carlo_cpp/tools/run_marseille_paper_batch.py
```

