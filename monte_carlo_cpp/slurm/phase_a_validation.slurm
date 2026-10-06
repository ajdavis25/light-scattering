#!/bin/bash
#SBATCH --job-name=phase_a_validation
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/phase_a_validation-%j.out
#
# Phase A (PAPER_READINESS_PLAN_2026-07-26.md): regenerate the validation-suite
# artifacts with the UNMODIFIED solver (build_cluster_gcc8, pre-fix code) so the
# historically claimed benchmark numbers are reproduced/reconciled, per audit
# finding M-2. Runs Zawada single smoke, Zawada multiple smoke, then the
# default gate, archiving validation_report__<label>.json per stage.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export RUNNER=monte_carlo_cpp/build_cluster_gcc8/ValidationRunner

echo "phase_a_slurm job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh
