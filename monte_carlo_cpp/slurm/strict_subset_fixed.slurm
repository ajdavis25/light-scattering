#!/bin/bash
#SBATCH --job-name=marseille_strict_subset_fixed
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/strict_subset_fixed-%j.out
#
# Phase B5 confirmation (PAPER_READINESS_PLAN_2026-07-26.md): strict-subset
# Marseille run with the FIXED build (build_fix): incoming-light polarization
# convention + unbiased mean higher-order estimator. Uses the resumable
# exact-direction batched runner with per-direction checkpointing so walltime
# kills lose at most one higher-order block. Offline predictions to compare
# against: raw median AoP error ~14 deg (was 43), masked shape RMSE ~0.12.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}

echo "strict_subset_fixed_slurm job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-16}"
"$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
    monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_strict_subset_fixed.cfg \
    --batch-size "${SLURM_CPUS_PER_TASK:-16}" \
    --higher-order-block-size 1 \
    --resume
