#!/bin/bash
#SBATCH --job-name=depol_lowguard_fix2
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/depol_lowguard_fix2-%j.out
#
# Depolarization-gap task #20 (2026-08-09), build_fix2 leg: low-order guard
# ladder (guard=3/4/6; needs the floor lowered 8 -> 2, build_fix2 only) plus
# the King-factor probe (rayleigh_depolarization_factor=0.0279; knob exists
# only in build_fix2 — older binaries silently ignore it). All on the
# Marseille 48-direction subset with chi-fix on, cap=1. build_fix2 verified:
# ctest 5/5, analytic DoLP(90)=(1-rho)/(1+rho) to 1e-16, rho=0 bit-identical.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export MONTE_CARLO_BUILD_DIR="$REPO_ROOT/monte_carlo_cpp/build_fix2"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
BASE=monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_branchcap_subset_chifix

echo "depol_lowguard_fix2 job=${SLURM_JOB_ID} node=$(hostname) workers=${SLURM_CPUS_PER_TASK:-32}"
for tag in guard3 guard4 guard6 depol0279; do
  echo "=== ${tag} start $(date -u +%H:%M:%SZ) ==="
  "$PYTHON_BIN" monte_carlo_cpp/tools/run_measurement_case_batched.py \
      "${BASE}_${tag}.cfg" \
      --batch-size "${SLURM_CPUS_PER_TASK:-32}" \
      --higher-order-block-size 1 \
      --resume
  echo "=== ${tag} done $(date -u +%H:%M:%SZ) ==="
done
echo "depol_lowguard_fix2 all done"
