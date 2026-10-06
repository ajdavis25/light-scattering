#!/bin/bash
#SBATCH --job-name=aero_loading_scan
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/aero_loading_scan-%j.out
#
# Task #21 (2026-08-10): aerosol-depolarization scan on the Marseille
# 48-direction corrected-chain subset. The bundled aerosol matrices are
# perfect-sphere Mie (f22 = f11 exactly), so aerosol events cannot
# depolarize; aerosol_depolarization_f22_ratio scales f22/f33/f44 at aerosol
# events (f11/f12 untouched -> intensity and single scatter unchanged).
# Scan 0.9/0.8/0.7/0.6 alone + the 0.7+King(0.0279) stacked candidate.
# Baselines: chifix subset +0.1938 signed median; King-only +0.1705.
# build_fix3: ctest 5/5; bit-identical to build_fix on the capped tiny case.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix3"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
BASE=monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_branchcap_subset_chifix

echo "aero_loading_scan job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
for tag in aerox2 aerox3 aerox3_full; do
  echo "=== ${tag} start $(date -u +%H:%M:%SZ) ==="
  "$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
      "${BASE}_${tag}.cfg" \
      --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
      --higher-order-block-size 1 \
      --resume
  echo "=== ${tag} done $(date -u +%H:%M:%SZ) ==="
done
echo "aero_loading_scan all done"
