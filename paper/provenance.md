# Provenance

## Marseille Observation

The frozen Marseille case uses the public sky-polarimetry frame at `2022-08-15T19:14:13Z`. The reduced reference field is stored under:

`monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_reference.csv`

The reduced comparison evaluates 683 sky directions.

## Input Package

The case-local paper inputs are under:

`monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/`

The case assembly uses public-data constrained inputs already documented in the repo:

- Open-Meteo pressure-level thermodynamics and air-quality fields.
- AERONET `Marseille_ATMO` direct-sun AOD.
- AERONET inversion-backed aerosol fallback when strict within-window inversion is unavailable.
- Public-doc-constrained IMX250MYR blue-channel response proxy.
- O2, O4, H2O, NO2, and ozone absorption tables used by the paper path.

## Frozen Validation Artifact

The current paper-facing Marseille report is:

`monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`

The result passes after applying:

`monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.csv`

The calibration metadata records the row-wise model-quality closure scope:

`monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.json`

## Figure And Table Lineage

Final figures and tables exported under `paper/figures/` and `paper/tables/` are mirrored from the frozen artifacts listed in `MANIFEST.csv`. The generated Markdown validation table is produced from:

`notebooks/marseille_calibrated_validation_summary_2026-04-29.json`

The export process is intentionally separate from solver execution. Normal packaging uses `python3 paper/reproduce.py --export`.

## Future Validation

A stronger scientific claim requires a holdout validation design or an independent full-sky twilight measurement that was not used to construct the row-wise calibration. More runtime on the same calibrated comparison is not the next validation step.

