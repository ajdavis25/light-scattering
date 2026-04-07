#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)
WORKER_SLURM_SCRIPT="${SCRIPT_DIR}/marseille_shard_worker_arc.slurm"
MERGE_SLURM_SCRIPT="${SCRIPT_DIR}/marseille_shard_merge_arc.slurm"
PREPARE_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/tools/prepare_measurement_case_shards.py"

PARTITION=${PARTITION:-compute1}
MERGE_PARTITION=${MERGE_PARTITION:-${PARTITION}}
JOB_NAME=${JOB_NAME:-marseille_strict_sharded}
TIME_LIMIT=${TIME_LIMIT:-18:00:00}
CPUS_PER_TASK=${CPUS_PER_TASK:-80}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-gcc/14.2.0}
MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG:-${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg}
CASE_ID_SUFFIX=${CASE_ID_SUFFIX:-arc_sharded}
SHARDS=${SHARDS:-8}
ARRAY_MAX=${ARRAY_MAX:-${SHARDS}}
HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-4}
BATCH_SIZE=${BATCH_SIZE:-${CPUS_PER_TASK}}
LAYOUT=${LAYOUT:-strided}

RESULTS_DIR="${REPO_ROOT}/monte_carlo_cpp/results"
SLURM_LOG_DIR="${RESULTS_DIR}/slurm"
REPORT_DIR="${RESULTS_DIR}/measurement_case_reports"
mkdir -p "${SLURM_LOG_DIR}" "${REPORT_DIR}"

BASE_CASE_ID=$(
    awk -F= '
        $1 ~ /^[[:space:]]*case_id[[:space:]]*$/ {
            value = $2
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            print value
            exit
        }
    ' "${MEASUREMENT_CONFIG}"
)
if [[ -z "${BASE_CASE_ID}" ]]; then
    echo "Unable to resolve case_id from ${MEASUREMENT_CONFIG}" >&2
    exit 1
fi

PARENT_CASE_ID="${BASE_CASE_ID}"
if [[ -n "${CASE_ID_SUFFIX}" ]]; then
    PARENT_CASE_ID="${BASE_CASE_ID}__${CASE_ID_SUFFIX}"
fi
MANIFEST_PATH="${REPORT_DIR}/_sharded/${PARENT_CASE_ID}/manifest.json"

"${PYTHON_BIN}" "${PREPARE_SCRIPT}" \
    "${MEASUREMENT_CONFIG}" \
    --shards "${SHARDS}" \
    --layout "${LAYOUT}" \
    --case-id-suffix "${CASE_ID_SUFFIX}" \
    --manifest "${MANIFEST_PATH}"

ACTUAL_SHARDS=$(
    "${PYTHON_BIN}" - <<'PY' "${MANIFEST_PATH}"
import json
import sys
from pathlib import Path
manifest = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
print(manifest["shard_count"])
PY
)
if (( ACTUAL_SHARDS < 1 )); then
    echo "No shards were created for ${MEASUREMENT_CONFIG}" >&2
    exit 1
fi
if (( ARRAY_MAX > ACTUAL_SHARDS )); then
    ARRAY_MAX=${ACTUAL_SHARDS}
fi

ARRAY_JOB_ID=$(
    sbatch --parsable \
        --partition="${PARTITION}" \
        --job-name="${JOB_NAME}" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="${CPUS_PER_TASK}" \
        --array="0-$((ACTUAL_SHARDS - 1))%${ARRAY_MAX}" \
        --time="${TIME_LIMIT}" \
        --output="${SLURM_LOG_DIR}/%x-%A_%a.out" \
        --error="${SLURM_LOG_DIR}/%x-%A_%a.err" \
        --export=ALL,REPO_ROOT="${REPO_ROOT}",MANIFEST_PATH="${MANIFEST_PATH}",PYTHON_BIN="${PYTHON_BIN}",GCC_MODULE="${GCC_MODULE}",BATCH_SIZE="${BATCH_SIZE}",HIGHER_ORDER_BLOCK_SIZE="${HIGHER_ORDER_BLOCK_SIZE}" \
        "${WORKER_SLURM_SCRIPT}"
)

MERGE_JOB_ID=$(
    sbatch --parsable \
        --dependency="afterok:${ARRAY_JOB_ID}" \
        --partition="${MERGE_PARTITION}" \
        --job-name="${JOB_NAME}_merge" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=1 \
        --time="01:00:00" \
        --output="${SLURM_LOG_DIR}/%x-%j.out" \
        --error="${SLURM_LOG_DIR}/%x-%j.err" \
        --export=ALL,REPO_ROOT="${REPO_ROOT}",MANIFEST_PATH="${MANIFEST_PATH}",PYTHON_BIN="${PYTHON_BIN}",GCC_MODULE="${GCC_MODULE}" \
        "${MERGE_SLURM_SCRIPT}"
)

ARRAY_JOB_ID_FILE="${SLURM_LOG_DIR}/${JOB_NAME}.latest_array_jobid"
MERGE_JOB_ID_FILE="${SLURM_LOG_DIR}/${JOB_NAME}.latest_merge_jobid"
printf '%s\n' "${ARRAY_JOB_ID}" > "${ARRAY_JOB_ID_FILE}"
printf '%s\n' "${MERGE_JOB_ID}" > "${MERGE_JOB_ID_FILE}"

echo "submitted_array_job_id=${ARRAY_JOB_ID}"
echo "submitted_merge_job_id=${MERGE_JOB_ID}"
echo "partition=${PARTITION}"
echo "merge_partition=${MERGE_PARTITION}"
echo "job_name=${JOB_NAME}"
echo "config=${MEASUREMENT_CONFIG}"
echo "manifest=${MANIFEST_PATH}"
echo "parent_case_id=${PARENT_CASE_ID}"
echo "requested_shards=${SHARDS}"
echo "actual_shards=${ACTUAL_SHARDS}"
echo "array_max_concurrency=${ARRAY_MAX}"
echo "cpus_per_task=${CPUS_PER_TASK}"
echo "batch_size=${BATCH_SIZE}"
echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
echo "array_job_id_file=${ARRAY_JOB_ID_FILE}"
echo "merge_job_id_file=${MERGE_JOB_ID_FILE}"
echo "merge_dependency=afterok:${ARRAY_JOB_ID}"
echo "parent_partial_rows_csv=${REPORT_DIR}/${PARENT_CASE_ID}_batched_partial_rows.csv"
echo "parent_progress_json=${REPORT_DIR}/${PARENT_CASE_ID}_batched_progress.json"
echo "parent_comparison_csv=${REPORT_DIR}/${PARENT_CASE_ID}_comparison.csv"
