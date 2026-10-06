#!/bin/bash
#SBATCH --job-name=rozenberg_discrim
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=04:00:00
#SBATCH --output=/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/rozenberg_discrim-%j.out
#
# Task #21 (2026-08-09): input discriminators for the twilight MS-bulk deficit.
# Three Rozenberg variants, each ~20 s: noo3 (ozone zeroed — bounds ozone
# suppression of the twilight field along tangent paths), noaero (aerosol
# zeroed — aerosol share of the diffuse sky), aerox4 (aerosol x4 —
# enhancement sensitivity). Note the model profile is aerosol-RICHER aloft
# than a true h^-6 law, so the aerosol-deficit hypothesis is already weak;
# refraction (shadow height 35 -> 29 km at SZA 96, ~2.2x illuminated
# molecular density) is the leading physics candidate and needs a code
# change, not a config run.
set -euo pipefail

REPO_ROOT=/work/vmo703/light-scattering
cd "$REPO_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-36}"
export RUNNER="$REPO_ROOT/monte_carlo_cpp/build_fix/ValidationRunner"

echo "rozenberg_discrim job=${SLURM_JOB_ID} node=$(hostname) omp=${OMP_NUM_THREADS}"
bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    split_rozenberg_noo3=monte_carlo_cpp/config/validation_split_rozenberg_noo3.cfg \
    split_rozenberg_noaero=monte_carlo_cpp/config/validation_split_rozenberg_noaero.cfg \
    split_rozenberg_aerox4=monte_carlo_cpp/config/validation_split_rozenberg_aerox4.cfg || true
echo "rozenberg_discrim done"
