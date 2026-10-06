#!/bin/bash
#SBATCH --job-name=validation_default_gate
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=9-00:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/validation_default_gate-%j.out
#
# The full unmodified default gate (convergence + DISORT + IPRT + Rozenberg +
# Koomen in ONE report) -- the artifact the historical claims came from. Job
# 778296 showed the convergence stage alone exceeds 48h at 16 threads (the sky
# loop IS OpenMP-parallel over 684 directions, so more threads scale): compute2
# allows 10-day walltime and this asks for 9. Rough core-hour extrapolation
# says ~6 days at 36 threads. stdbuf in the driver line-buffers the runner
# stdout, so the two convergence sky runs are visible as two bursts of
# "Evaluating zenith row" lines in default_clear_sky_stdout.log.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER=monte_carlo_cpp/build_cluster_gcc8/ValidationRunner

echo "validation_default_gate job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    default_clear_sky=monte_carlo_cpp/config/default_clear_sky.cfg
