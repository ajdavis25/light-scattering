#!/bin/bash
#SBATCH --job-name=depol_guard_ladder
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/depol_guard_ladder-%j.out
#
# Depolarization-gap task #20 (2026-08-09): max_events_guard ladder on the
# Marseille 48-direction subset with the corrected chain (chi-fix on, cap=1).
# Measures total-sky DoLP as a function of maximum scattering order to locate
# where the over-polarization (+0.194 signed median at guard=64, job 781658)
# saturates. guard=64 baseline already exists (781658); this job adds
# guard=8/16/32 on the same build_fix binary. The runtime floors the guard at
# 8, so orders below 8 need the build_fix2 variant (separate job).
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}

echo "depol_guard_ladder job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
for g in 8 16 32; do
  echo "=== guard ${g} start $(date -u +%H:%M:%SZ) ==="
  "$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
      "monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_branchcap_subset_chifix_guard${g}.cfg" \
      --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
      --higher-order-block-size 1 \
      --resume
  echo "=== guard ${g} done $(date -u +%H:%M:%SZ) ==="
done
echo "depol_guard_ladder all done"
