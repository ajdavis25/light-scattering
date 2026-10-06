#!/bin/bash
#SBATCH --job-name=validation_default_capped
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=2-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/validation_default_capped-%j.out
#
# Full default gate under the declared cap=1 deviation
# (config/validation_default_capped.cfg). Submitted 2026-08-03 while 779982
# (the unmodified gate) sits silent 4.8 d into its 9 d wall: the uncapped
# SZA-97 all-orders convergence sky is the runtime-confirmed non-terminating
# recursion regime, so 779982 is expected to TIMEOUT empty; it is left running
# as the legacy control. Scaling from 779983 (140 s/direction single-thread at
# 25 bands / 256 photons, cap=1), this capped gate's two convergence skies
# (684 dirs x 46 bands, 128+256 photons) should take ~2 h at 36 threads, plus
# the four case stages. Convergence metrics land in
# results/validation/validation_report__default_clear_sky_capped.json.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER=monte_carlo_cpp/build_cluster_gcc8/ValidationRunner

echo "validation_default_capped job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    default_clear_sky_capped=monte_carlo_cpp/config/validation_default_capped.cfg
