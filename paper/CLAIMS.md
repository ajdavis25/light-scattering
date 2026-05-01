# Claim Discipline

## Allowed Core Claim

The repo provides a reproducible polarized twilight radiative-transfer pipeline with benchmark checks and a calibrated full-field Marseille validation artifact. The Marseille full-field comparison passes as calibrated pipeline validation after applying the frozen row-wise measurement-model calibration.

## Required Limitation

The Marseille result is calibrated pipeline validation. It should not be interpreted as independent raw first-principles validation. The row-wise measurement-model calibration was derived on the same Marseille comparison grid, so a stronger physical claim requires a holdout sky measurement, a calibration split, or an independent public full-sky twilight measurement.

## Allowed Supporting Claims

- The C++ solver is the scientific source of truth for the current pipeline.
- The Python layer is downstream plotting and comparison support.
- The frozen Marseille package validates data reduction, comparison-grid plumbing, metric extraction, plotting, and calibrated measurement-path closure.
- The exported figures and tables are derived from frozen artifacts listed in `MANIFEST.csv`.
- Historical raw Marseille mismatches remain useful as diagnostic and limitation context.

## Forbidden Or Unqualified Claims

- Do not claim independent raw first-principles closure for Marseille.
- Do not claim raw predictive validation of the Marseille twilight polarization field.
- Do not claim that the calibrated Marseille pass proves first-principles below-horizon twilight polarization closure.
- Do not present order-fraction diagnostic plots as proof of raw model closure.
- Do not imply that another identical expensive Marseille rerun would strengthen the calibrated validation claim by itself.

## Referee-Response Wording

We separate the Marseille result into calibrated pipeline validation and raw physics validation. The reported full-field Marseille artifact applies a frozen row-wise measurement-model calibration and passes the configured gates over 683 sky directions. This demonstrates that the data reduction, comparison grid, metric extraction, plotting, and calibrated measurement pathway are internally closed. We do not present this as independent first-principles twilight closure; that stronger claim would require a holdout sky measurement or a calibration split not used to construct the row-wise mapping.

