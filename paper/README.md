# Marseille Calibrated Pipeline Validation Paper Package

This directory is the publication-facing package for the Marseille calibrated pipeline validation paper.

The defensible paper claim is:

> The repo provides a reproducible polarized twilight radiative-transfer pipeline with benchmark checks and a calibrated full-field Marseille validation artifact. The Marseille full-field comparison passes as calibrated pipeline validation after applying the frozen row-wise measurement-model calibration.

The Marseille result must not be described as independent raw first-principles closure. The pass validates the calibrated data path, comparison grid, metric extraction, plotting, and configured gate passage after frozen row-wise calibration.

## Start Here

- `CLAIMS.md`: allowed claims, forbidden claims, limitation wording, and referee-response language.
- `MANIFEST.csv`: frozen source artifacts, export targets, and SHA256 checksums.
- `provenance.md`: source data, generated artifacts, and validation lineage.
- `text_snippets.md`: manuscript-ready methods, results, limitation, table, and caption text.
- `reproduce.py`: verification and export CLI.

## Commands

Verify the frozen package without running the solver:

```bash
python3 paper/reproduce.py --verify-only
```

Export final paper-facing figures and tables from frozen artifacts:

```bash
python3 paper/reproduce.py --export
```

Regenerate Marseille plots from the frozen report/comparison CSV, then export:

```bash
python3 paper/reproduce.py --regenerate-plots
```

Plot regeneration requires an interpreter new enough to import the repo plotting modules plus `numpy` and `matplotlib`. If those dependencies are not available, use `--export`; it mirrors the already-frozen figures without touching the solver.

The normal paper workflow does not rerun the expensive Marseille solver. Strict Marseille reruns are out of scope unless the calibration, raw model, input package, or paper claim changes.

## Canonical Source Of Truth

The canonical frozen Marseille artifacts stay in their existing repo locations under `notebooks/`, `monte_carlo_cpp/results/measurement_case_reports/`, and `monte_carlo_cpp/data/paper_cases/`. This `paper/` directory mirrors the final publication-facing figures and tables, and records the exact source artifacts used for those exports.
