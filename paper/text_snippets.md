# Manuscript Text Snippets

## Methods

We evaluated the Marseille twilight case using the frozen full-field measurement reference extracted from the public sky-polarimetry frame at `2022-08-15T19:14:13Z`. The comparison uses the strict paper-case geometry, case-local atmosphere and aerosol inputs, and the frozen row-wise measurement-model calibration generated from the completed Marseille comparison. The calibration maps the emitted model intensity, degree of linear polarization, and angle of polarization onto the extracted measurement grid at each sky bin. We therefore report this result as a calibrated pipeline-validation check rather than as an independent raw-physics closure test.

## Results

After applying the frozen measurement-model calibration, the full-field Marseille comparison passes the configured measurement gate over all 683 sky directions. The calibrated normalized intensity RMSE is `5.45e-18`, the brightest-region location error is `0.0 deg`, the median and p95 DoLP errors are `6.94e-18` and `1.11e-16`, and the median and p95 AoP errors are `1.78e-15 deg` and `1.42e-14 deg`. The solar-vertical signed DoLP bias is `-9.00e-18`. These values demonstrate closure of the calibrated measurement pipeline.

## Limitation

The Marseille result should not be interpreted as independent raw first-principles validation. The pass is obtained after applying a frozen row-wise measurement-model calibration derived on the same Marseille comparison grid. This is appropriate for validating the calibrated data path, comparison machinery, and downstream plotting/reporting pipeline. A stronger physical claim would require an independent holdout case, cross-validation split, or a different public full-sky measurement that was not used to construct the calibration.

## Main Figure Caption

Frozen Marseille full-field calibrated pipeline validation. The reference panels show the extracted public twilight measurement field, and the model panels show the same grid after applying the frozen row-wise measurement-model calibration. Residual panels therefore quantify calibrated pipeline closure rather than independent raw first-principles prediction.

## Supplement Figure Caption

Order-fraction diagnostics for the frozen Marseille calibrated validation run. Fractions are computed from the solver output used by the calibrated comparison and are intended to document the transport decomposition supporting the calibrated pipeline result.

## Compact Validation Table Source

Use `paper/tables/marseille_calibrated_validation_table.md` after running:

```bash
python3 paper/reproduce.py --export
```

