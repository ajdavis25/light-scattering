#!/bin/bash
#SBATCH --job-name=marseille_subset_aerox4
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/branchcap_subset_chifix_aerox4-%j.out
#
# Marseille 48-direction branchcap subset with event_frame_chi_sign_fix=true.
# Identical to the 779983 run (which completed 48/48 in 209 s) except for the
# chi-sign fix knob. Purpose: test whether the eventMuellerMatrix chi-in
# handedness bug (found 2026-08-01; fixing it makes ALL FOUR IPRT A1 vector
# metrics reproduce their historical claims almost digit-for-digit) is also
# the mechanism behind the persistent Marseille DoLP over-polarization.
#
# Key readouts vs the 779983 baseline (case ..._branchcap_subset_fixed):
#   signed DoLP err median: baseline +0.1431 (mean +0.1303). If the chi bug is
#     the mechanism, this should drop markedly, especially in the
#     solar-vertical midzen region (baseline bias +0.307, second_rr_frac 0.56
#     dominant there).
#   raw median |AoP err|: baseline 22.29 deg -- should not regress materially.
#   masked shape RMSE: baseline 0.2085 -- direction of change unconstrained.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}

echo "branchcap_subset_chifix_slurm job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
"$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
    monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_branchcap_subset_chifix_aerox4.cfg \
    --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
    --higher-order-block-size 1 \
    --resume
