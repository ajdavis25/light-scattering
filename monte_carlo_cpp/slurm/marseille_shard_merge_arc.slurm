#!/bin/bash
set -euo pipefail

if [[ -f /etc/profile.d/lmod.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/lmod.sh
elif [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

REPO_ROOT=${REPO_ROOT:?REPO_ROOT must be exported by the submit helper}
MANIFEST_PATH=${MANIFEST_PATH:?MANIFEST_PATH must be exported by the submit helper}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-gcc/14.2.0}

if ! command -v module >/dev/null 2>&1; then
    echo "module command is unavailable in this Slurm shell" >&2
    exit 1
fi

module purge
module load "${GCC_MODULE}"

export PYTHONUNBUFFERED=1

cd "${REPO_ROOT}"

MERGE_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/tools/merge_measurement_case_shards.py"
for required_path in "${PYTHON_BIN}" "${MERGE_SCRIPT}" "${MANIFEST_PATH}"; do
    if [[ ! -e "${required_path}" ]]; then
        echo "Required path is missing: ${required_path}" >&2
        exit 1
    fi
done

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "host=$(hostname)"
echo "slurm_job_id=${SLURM_JOB_ID:-}"
echo "manifest_path=${MANIFEST_PATH}"

"${PYTHON_BIN}" -u "${MERGE_SCRIPT}" "${MANIFEST_PATH}"

echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
