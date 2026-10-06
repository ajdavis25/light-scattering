#!/bin/bash
#SBATCH --job-name=phase_a_split_cases
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=3-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/phase_a_split_cases-%j.out
#
# Phase A continuation after job 778296 timed out inside the convergence stage:
# the default gate runs convergence FIRST (2 full-sky runs of the top-level cfg),
# so none of the four benchmark/measurement case artifacts ever landed. These
# wrappers shrink the top-level sky to 1x1 bins (the case configs drive their
# own runs, proven by the zawada smokes) so each case artifact lands on its own.
# Sequential fastest-first with per-stage archiving; koomen runs in a sibling job.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER=monte_carlo_cpp/build_cluster_gcc8/ValidationRunner

echo "phase_a_split_cases job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    split_disort=monte_carlo_cpp/config/validation_split_disort.cfg \
    split_iprt=monte_carlo_cpp/config/validation_split_iprt.cfg \
    split_rozenberg=monte_carlo_cpp/config/validation_split_rozenberg.cfg
