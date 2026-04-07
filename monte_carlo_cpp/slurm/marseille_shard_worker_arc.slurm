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
BATCH_SIZE=${BATCH_SIZE:-${SLURM_CPUS_PER_TASK:-1}}
HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-4}
SHARD_INDEX=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

if ! command -v module >/dev/null 2>&1; then
    echo "module command is unavailable in this Slurm shell" >&2
    exit 1
fi

module purge
module load "${GCC_MODULE}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

cd "${REPO_ROOT}"

RUNNER_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/tools/run_measurement_case_batched.py"
for required_path in "${PYTHON_BIN}" "${RUNNER_SCRIPT}" "${MANIFEST_PATH}"; do
    if [[ ! -e "${required_path}" ]]; then
        echo "Required path is missing: ${required_path}" >&2
        exit 1
    fi
done

readarray -t shard_meta < <(
    "${PYTHON_BIN}" - <<'PY' "${MANIFEST_PATH}" "${SHARD_INDEX}"
import json
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
shard = manifest["shards"][int(sys.argv[2])]
print(shard["case_id"])
print(shard["config"])
print(shard["direction_count"])
print(shard["partial_rows_csv"])
print(shard["progress_json"])
PY
)

CASE_ID=${shard_meta[0]}
SHARD_CONFIG=${shard_meta[1]}
DIRECTION_COUNT=${shard_meta[2]}
PARTIAL_ROWS=${shard_meta[3]}
PROGRESS_JSON=${shard_meta[4]}

LOCAL_BATCH_SIZE=${BATCH_SIZE}
if (( DIRECTION_COUNT < LOCAL_BATCH_SIZE )); then
    LOCAL_BATCH_SIZE=${DIRECTION_COUNT}
fi
if (( LOCAL_BATCH_SIZE < 1 )); then
    echo "No directions assigned to shard ${SHARD_INDEX}; skipping."
    exit 0
fi

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "host=$(hostname)"
echo "slurm_job_id=${SLURM_JOB_ID:-}"
echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-}"
echo "case_id=${CASE_ID}"
echo "manifest_path=${MANIFEST_PATH}"
echo "config=${SHARD_CONFIG}"
echo "direction_count=${DIRECTION_COUNT}"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-}"
echo "batch_size=${LOCAL_BATCH_SIZE}"
echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
echo "partial_rows_csv=${PARTIAL_ROWS}"
echo "progress_json=${PROGRESS_JSON}"

srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}" --cpu-bind=cores \
    "${PYTHON_BIN}" -u "${RUNNER_SCRIPT}" \
    "${SHARD_CONFIG}" \
    --batch-size "${LOCAL_BATCH_SIZE}" \
    --higher-order-block-size "${HIGHER_ORDER_BLOCK_SIZE}" \
    --resume

echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
