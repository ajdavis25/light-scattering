#!/bin/bash
#SBATCH --job-name=marseille_branchcap_subset_fixed
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=6-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/branchcap_subset_fixed-%j.out
#
# Phase B5 confirmation, second attempt. Job 778297 (strict_subset tier) proved
# computationally infeasible: that tier REMOVES higher_order_recursive_branch_cap,
# so 3-way twilight branching explodes and 0/48 directions finished in 48h.
# This run uses the SAME solver budget as the frozen __full_branchcap_robust_r2
# production run (256 photons/bin, branch cap 1, robust groups 16) -- the tier
# every offline prediction was computed from -- restricted to the 48-direction
# subset, with the two config-gated fixes enabled.
#
# Subset-specific offline predictions (recomputed 2026-07-29 within the 48-dir
# field; full-field numbers do NOT transfer):
#   raw median |AoP err|: 52.5 deg -> 21.2 deg with incoming_light convention
#   masked shape RMSE: ~0.253 -> ~0.255 (estimator swap does NOT improve shape
#     on this subset -- do not use shape RMSE as the estimator pass criterion;
#     use per-direction consistency with frozen totals re-aggregated as mean)
#   DoLP signed err median: ~ +0.134 persists (genuine physics gap)
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}

echo "branchcap_subset_fixed_slurm job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
"$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
    monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_branchcap_subset_fixed.cfg \
    --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
    --higher-order-block-size 1 \
    --resume
