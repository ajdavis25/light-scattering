#!/bin/bash
#SBATCH --job-name=validation_capped_chifix
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=04:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/validation_capped_chifix-%j.out
#
# Chi-fix adoption prep, task #23 (2026-08-09): the full default gate under the
# declared cap=1 deviation WITH event_frame_chi_sign_fix=true — the candidate
# production gate. Baseline: the legacy capped gate (783441) ran everything in
# 21 s; expect similar here. Uses build_fix (knows the chi-fix knob; the Apr-27
# build_cluster_gcc8 binary would silently ignore it). Report lands as
# validation_report__default_clear_sky_capped_chifix.json.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER="$REPO_ROOT/monte_carlo_cpp/build_fix/ValidationRunner"

echo "validation_capped_chifix job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    default_clear_sky_capped_chifix=monte_carlo_cpp/config/validation_default_capped_chifix.cfg
