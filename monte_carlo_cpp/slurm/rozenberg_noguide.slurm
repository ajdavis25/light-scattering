#!/bin/bash
#SBATCH --job-name=rozenberg_noguide
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=1-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/rozenberg_noguide-%j.out
#
# Task #21 (2026-08-09): Rozenberg MS-deficit discriminator. The fresh split
# rerun shows the model sky falling off the horizon peak 3-15x faster than the
# Rozenberg reference, single-scatter dominated where the real twilight sky is
# multiple-scatter dominated. This run disables the twilight higher-order
# guided continuation (plain single-branch phase sampling, the vanilla
# unbiased estimator) to test whether the guided cap=1 estimator is biased low
# in its multiple-scatter yield. Compare higher_frac and normalized profile
# against validation_report__split_rozenberg.json / its comparison CSV.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER="$REPO_ROOT/monte_carlo_cpp/build_fix/ValidationRunner"

echo "rozenberg_noguide job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    split_rozenberg_noguide=monte_carlo_cpp/config/validation_split_rozenberg_noguide.cfg
