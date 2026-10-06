#!/bin/bash
#SBATCH --job-name=gate_v2
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=04:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/gate_v2-%j.out
#
# V2 production gate, task #23 (2026-08-10): capped chi-fix gate with
# (a) hardened convergence metric (convergence_low_order_metric=true — peaks
#     from first+second order; heavy-tailed higher-order excluded; flux total),
# (b) h-6-profile measurement cases with King 0.0279 (Rozenberg + Koomen both
#     get the true h-6 profile; same input archaeology, applied coherently),
# (c) benchmarks unchanged (published spec, no King).
# Uses the PINNED build_v2freeze binary (never rebuilt during runs; the
# build_fix3 in-place rebuild killed 788177 with ETXTBSY).
# Report lands as validation_report__default_clear_sky_capped_v2.json.
# Expected wall: ~21 s class (783441 baseline). Slurm FAILED state is the
# expected nonzero exit if overall_pass=false (Koomen DoLP declared fail).
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER="$REPO_ROOT/monte_carlo_cpp/build_v2freeze/ValidationRunner"

echo "gate_v2 job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    default_clear_sky_capped_v2=monte_carlo_cpp/config/validation_default_capped_v2.cfg
