#!/bin/bash
#SBATCH --job-name=phase_a_split_koomen
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=5-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/phase_a_split_koomen-%j.out
#
# Koomen is the heaviest of the four default-gate cases (19x36 bins x 1024
# photons, all orders, sun at h=-6 deg), so it gets its own walltime instead of
# queueing behind the other three in phase_a_split_cases.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER=monte_carlo_cpp/build_cluster_gcc8/ValidationRunner

echo "phase_a_split_koomen job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    split_koomen=monte_carlo_cpp/config/validation_split_koomen.cfg
