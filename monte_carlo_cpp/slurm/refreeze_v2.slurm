#!/bin/bash
#SBATCH --job-name=refreeze_v2
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/refreeze_v2-%j.out
#
# Task #23 (2026-08-10): v2 physics-only re-freeze runs.
# Stack: chi-fix + King 0.0279 + AERONET at-site day profile (BL x1.3796,
# column AOD550 = 0.155 from Marseille_ATMO, the station bundled in the frozen
# case dir; L2.0-confirmed). f22 excluded (sensitivity only).
# Order: 48-dir subset sanity (~4 min), then 683-dir full field (~55 min).
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

# Pinned freeze binary: NEVER rebuilt during runs (build_fix3 rebuild mid-job
# caused ETXTBSY on 788177 — dev builds and frozen-tier binaries must differ).
export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_v2freeze"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
BASE=monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement

echo "refreeze_v2 job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
for tag in branchcap_subset_v2 full_v2; do
  echo "=== ${tag} start $(date -u +%H:%M:%SZ) ==="
  "$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
      "${BASE}_${tag}.cfg" \
      --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
      --higher-order-block-size 1 \
      --resume
  echo "=== ${tag} done $(date -u +%H:%M:%SZ) ==="
done
echo "refreeze_v2 all done"
