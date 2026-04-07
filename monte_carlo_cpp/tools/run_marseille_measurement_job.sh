#!/bin/bash
set -euo pipefail

if [[ -f /etc/profile.d/lmod.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/lmod.sh
elif [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

REPO_ROOT=${REPO_ROOT:?REPO_ROOT must be set by the scheduler submit wrapper}
MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG:?MEASUREMENT_CONFIG must be set by the scheduler submit wrapper}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-}
MONTE_CARLO_BUILD_DIR=${MONTE_CARLO_BUILD_DIR:-${REPO_ROOT}/monte_carlo_cpp/build_cluster_gcc8}
BATCH_SIZE=${BATCH_SIZE:-40}
HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-1}
SCHEDULER_FAMILY=${SCHEDULER_FAMILY:-unknown}

MODULE_STATUS=not_requested
if [[ -n "${GCC_MODULE}" ]]; then
    if command -v module >/dev/null 2>&1; then
        if module spider "${GCC_MODULE}" >/dev/null 2>&1; then
            module purge
            module load "${GCC_MODULE}"
            MODULE_STATUS="loaded:${GCC_MODULE}"
        else
            MODULE_STATUS="missing:${GCC_MODULE}"
        fi
    else
        MODULE_STATUS="module_command_unavailable:${GCC_MODULE}"
    fi
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export MONTE_CARLO_BUILD_DIR

ulimit -n 4096 || true

cd "${REPO_ROOT}"

RUNNER_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/tools/run_measurement_case_batched.py"
for required_path in "${PYTHON_BIN}" "${RUNNER_SCRIPT}" "${MEASUREMENT_CONFIG}"; do
    if [[ ! -e "${required_path}" ]]; then
        echo "Required path is missing: ${required_path}" >&2
        exit 1
    fi
done

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

if [[ -z "${CASE_ID}" ]]; then
    echo "Unable to resolve case_id from ${MEASUREMENT_CONFIG}" >&2
    exit 1
fi

RESULTS_DIR="${REPO_ROOT}/monte_carlo_cpp/results"
REPORT_DIR="${RESULTS_DIR}/measurement_case_reports"
RUN_LOG="${REPORT_DIR}/${CASE_ID}_batched.log"
PROGRESS_JSON="${REPORT_DIR}/${CASE_ID}_batched_progress.json"
PARTIAL_ROWS="${REPORT_DIR}/${CASE_ID}_batched_partial_rows.csv"

mkdir -p "${REPORT_DIR}"

on_exit() {
    local exit_code=$?
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "exit_code=${exit_code}"
}

trap on_exit EXIT

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "host=$(hostname)"
echo "scheduler_family=${SCHEDULER_FAMILY}"
echo "scheduler_job_id=${SLURM_JOB_ID:-${PBS_JOBID:-}}"
echo "scheduler_queue=${SLURM_JOB_PARTITION:-${PBS_QUEUE:-}}"
echo "cpus_requested=${SLURM_CPUS_PER_TASK:-${PBS_NP:-}}"
echo "repo_root=${REPO_ROOT}"
echo "measurement_config=${MEASUREMENT_CONFIG}"
echo "python_bin=${PYTHON_BIN}"
echo "module_status=${MODULE_STATUS}"
echo "monte_carlo_build_dir=${MONTE_CARLO_BUILD_DIR}"
echo "batch_size=${BATCH_SIZE}"
echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
echo "resume_mode=enabled"
echo "measurement_log=${RUN_LOG}"
echo "progress_json=${PROGRESS_JSON}"
echo "partial_rows_csv=${PARTIAL_ROWS}"

"${PYTHON_BIN}" -u "${RUNNER_SCRIPT}" \
    "${MEASUREMENT_CONFIG}" \
    --batch-size "${BATCH_SIZE}" \
    --higher-order-block-size "${HIGHER_ORDER_BLOCK_SIZE}" \
    --resume
