#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)
SLURM_SCRIPT="${SCRIPT_DIR}/marseille_measurement_arc.slurm"

PARTITION=${PARTITION:-anantuabhg}
JOB_NAME=${JOB_NAME:-marseille_strict_measurement}
TIME_LIMIT=${TIME_LIMIT:-2-00:00:00}
CPUS_PER_TASK=${CPUS_PER_TASK:-80}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-gcc/14.2.0}
MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG:-${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg}
HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-32}
BATCH_SIZE=${BATCH_SIZE:-${CPUS_PER_TASK}}

RESULTS_DIR="${REPO_ROOT}/monte_carlo_cpp/results"
SLURM_LOG_DIR="${RESULTS_DIR}/slurm"
REPORT_DIR="${RESULTS_DIR}/measurement_case_reports"

mkdir -p "${SLURM_LOG_DIR}" "${REPORT_DIR}"

JOB_ID=$(
    sbatch --parsable \
        --partition="${PARTITION}" \
        --job-name="${JOB_NAME}" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="${CPUS_PER_TASK}" \
        --time="${TIME_LIMIT}" \
        --output="${SLURM_LOG_DIR}/%x-%j.out" \
        --error="${SLURM_LOG_DIR}/%x-%j.err" \
        --export=ALL,REPO_ROOT="${REPO_ROOT}",MEASUREMENT_CONFIG="${MEASUREMENT_CONFIG}",PYTHON_BIN="${PYTHON_BIN}",GCC_MODULE="${GCC_MODULE}",BATCH_SIZE="${BATCH_SIZE}",HIGHER_ORDER_BLOCK_SIZE="${HIGHER_ORDER_BLOCK_SIZE}" \
        "${SLURM_SCRIPT}"
)

JOB_ID_FILE="${SLURM_LOG_DIR}/${JOB_NAME}.latest_jobid"
printf '%s\n' "${JOB_ID}" > "${JOB_ID_FILE}"

CASE_ID=$(
    awk -F= '
        $1 ~ /^[[:space:]]*case_id[[:space:]]*$/ {
            value = $2
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            print value
            exit
        }
    ' "${MEASUREMENT_CONFIG}"
)

RUN_LOG="${REPORT_DIR}/${CASE_ID}_batched.log"
PROGRESS_JSON="${REPORT_DIR}/${CASE_ID}_batched_progress.json"
PARTIAL_ROWS="${REPORT_DIR}/${CASE_ID}_batched_partial_rows.csv"

echo "submitted_job_id=${JOB_ID}"
echo "partition=${PARTITION}"
echo "job_name=${JOB_NAME}"
echo "config=${MEASUREMENT_CONFIG}"
echo "cpus_per_task=${CPUS_PER_TASK}"
echo "batch_size=${BATCH_SIZE}"
echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
echo "slurm_stdout=${SLURM_LOG_DIR}/${JOB_NAME}-${JOB_ID}.out"
echo "slurm_stderr=${SLURM_LOG_DIR}/${JOB_NAME}-${JOB_ID}.err"
echo "measurement_log=${RUN_LOG}"
echo "progress_json=${PROGRESS_JSON}"
echo "partial_rows_csv=${PARTIAL_ROWS}"
echo "job_id_file=${JOB_ID_FILE}"
