#!/bin/bash
set -euo pipefail

if [[ -f /etc/profile.d/lmod.sh ]]; then
    # Ensure `module` is defined in non-login Slurm shells.
    # shellcheck disable=SC1091
    source /etc/profile.d/lmod.sh
elif [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

REPO_ROOT=${REPO_ROOT:?REPO_ROOT must be exported by the submit helper}
MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG:?MEASUREMENT_CONFIG must be exported by the submit helper}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-gcc/14.2.0}
BATCH_SIZE=${BATCH_SIZE:-${SLURM_CPUS_PER_TASK:-1}}
HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-32}

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

ulimit -n 4096 || true

cd "${REPO_ROOT}"

BUILD_DIR="${REPO_ROOT}/monte_carlo_cpp/build_current"
MEASUREMENT_RUNNER="${BUILD_DIR}/MeasurementCaseRunner"
VALIDATION_RUNNER="${BUILD_DIR}/ValidationRunner"
RUNNER_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/tools/run_measurement_case_batched.py"

for required_path in "${PYTHON_BIN}" "${MEASUREMENT_RUNNER}" "${VALIDATION_RUNNER}" "${RUNNER_SCRIPT}" "${MEASUREMENT_CONFIG}"; do
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
PROGRESS_JSON="${REPORT_DIR}/${CASE_ID}_batched_progress.json"
RUN_LOG="${REPORT_DIR}/${CASE_ID}_batched.log"
PARTIAL_ROWS="${REPORT_DIR}/${CASE_ID}_batched_partial_rows.csv"

mkdir -p "${REPORT_DIR}" "${RESULTS_DIR}/slurm"

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "host=$(hostname)"
echo "slurm_job_id=${SLURM_JOB_ID:-}"
echo "partition=${SLURM_JOB_PARTITION:-}"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-}"
echo "repo_root=${REPO_ROOT}"
echo "measurement_config=${MEASUREMENT_CONFIG}"
echo "measurement_runner=${MEASUREMENT_RUNNER}"
echo "validation_runner=${VALIDATION_RUNNER}"
echo "python_bin=${PYTHON_BIN}"
echo "gcc_module=${GCC_MODULE}"
echo "batch_size=${BATCH_SIZE}"
echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
echo "run_log=${RUN_LOG}"
echo "progress_json=${PROGRESS_JSON}"
echo "partial_rows_csv=${PARTIAL_ROWS}"

srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}" --cpu-bind=cores \
    "${PYTHON_BIN}" -u "${RUNNER_SCRIPT}" \
    "${MEASUREMENT_CONFIG}" \
    --batch-size "${BATCH_SIZE}" \
    --higher-order-block-size "${HIGHER_ORDER_BLOCK_SIZE}" \
    --resume

echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
