# light-scattering

## Project Description

light scattering code

## Current Marseille Validation Status

The current full-field Marseille pipeline-valid result is the calibrated frozen report at `monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`.

Use `notebooks/MARSEILLE_CALIBRATED_VALIDATION_STATUS_2026-04-28.md` as the downstream project note and `notebooks/MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md` for manuscript-ready wording, figure inventory, and validation-summary pointers. The Marseille gate passes after applying the frozen row-wise measurement-model calibration, and the result should not be described as independent raw first-principles closure.


## SMOOTHIE BINGUS
## SUIT BINGUS
## KINGUS


## Project Structure

light-scattering/
├── data/
│   ├── 
│   └── 
├── notebooks/
│   ├── plots/
│   ├── older_code_versions.ipynb
│   ├── tests_imports.ipynb
│   ├── test_scripts.ipynb
│   └── tests.ipynb
├── plots/
├── src/
│   ├── __init__.py
│   ├── intensity_vs_theta.py
│   ├── logging_setup.py
│   ├── main.py
│   ├── plane_parallel.py
│   ├── polarization_vs_theta.py
│   ├── single_scatter_plot.py
│   ├── utils.py
│   └── vector_operations.py
├── tests/
│   ├── __init__.py
│   ├── intensity_vs_theta_test.py
│   ├── plane_parallel_test.py
│   ├── polarization_vs_theta_test.py
│   ├── single_scatter_plot_test.py
│   └── vector_operations_test.py
├── __init__.py
└── README.md


## Unittests

To run unittests, open cmd and navigate to the root directory folder:
- Set your Python path by running an command similar to `set PYTHONPATH=c:/Users/ashton/Desktop/light-scattering/src`
- Run the command `pytest tests/whatever_script.py` (replace whatever_script.py to the desired .py to test)
