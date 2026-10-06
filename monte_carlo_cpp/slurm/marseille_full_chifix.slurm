#!/bin/bash
#SBATCH --job-name=marseille_full_chifix
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_full_chifix-%j.out
#
# Chi-fix adoption prep (task #23) + depolarization gap (task #20), 2026-08-09:
# full-field 683-direction Marseille run with the corrected chain (chi-fix on,
# incoming_light convention, mean estimator, cap=1). Gives (a) the
# adoption-candidate full-field metrics and (b) region-resolved signed-DoLP
# bias over the whole sky for the depolarization analysis. Scaling from the
# 48-direction subset (781658: 48/48 in 3m31s at 32 workers), expect ~50-60 min.
# Does NOT touch any frozen artifact (new case_id _full_chifix).
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}

echo "marseille_full_chifix job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
"$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
    monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_full_chifix.cfg \
    --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
    --higher-order-block-size 1 \
    --resume
echo "marseille_full_chifix done"
