#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)
SLURM_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/slurm/marseille_measurement_longhaul.slurm"
PBS_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/pbs/marseille_measurement_longhaul.pbs"

SCHEDULER=${SCHEDULER:-auto}
PARTITION=${PARTITION:-compute1}
QUEUE=${QUEUE:-${PARTITION}}
JOB_NAME=${JOB_NAME:-marseille_strict_measurement_longhaul}
TIME_LIMIT=${TIME_LIMIT:-2-00:00:00}
CPUS_PER_TASK=${CPUS_PER_TASK:-40}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-}
MONTE_CARLO_BUILD_DIR=${MONTE_CARLO_BUILD_DIR:-${REPO_ROOT}/monte_carlo_cpp/build_cluster_gcc8}
MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG:-${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg}
BATCH_SIZE=${BATCH_SIZE:-${CPUS_PER_TASK}}
HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-1}
SLURM_HINT=${SLURM_HINT:-nomultithread}

RESULTS_DIR="${REPO_ROOT}/monte_carlo_cpp/results"
REPORT_DIR="${RESULTS_DIR}/measurement_case_reports"
SLURM_LOG_DIR="${RESULTS_DIR}/slurm"
PBS_LOG_DIR="${RESULTS_DIR}/pbs"
mkdir -p "${REPORT_DIR}" "${SLURM_LOG_DIR}" "${PBS_LOG_DIR}"

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

RUN_LOG="${REPORT_DIR}/${CASE_ID}_batched.log"
PROGRESS_JSON="${REPORT_DIR}/${CASE_ID}_batched_progress.json"
PARTIAL_ROWS="${REPORT_DIR}/${CASE_ID}_batched_partial_rows.csv"

detect_scheduler() {
    if [[ "${SCHEDULER}" == "slurm" || "${SCHEDULER}" == "pbs" ]]; then
        printf '%s\n' "${SCHEDULER}"
        return
    fi
    if command -v sbatch >/dev/null 2>&1; then
        printf 'slurm\n'
        return
    fi
    if command -v qsub >/dev/null 2>&1; then
        printf 'pbs\n'
        return
    fi
    echo "No supported scheduler command found. Expected sbatch or qsub." >&2
    exit 1
}

slurm_time_to_pbs() {
    local slurm_time=$1
    if [[ "${slurm_time}" =~ ^([0-9]+)-([0-9]{2}):([0-9]{2}):([0-9]{2})$ ]]; then
        local days=${BASH_REMATCH[1]}
        local hours=${BASH_REMATCH[2]}
        local minutes=${BASH_REMATCH[3]}
        local seconds=${BASH_REMATCH[4]}
        printf '%d:%02d:%02d\n' \
            "$((10#${days} * 24 + 10#${hours}))" \
            "$((10#${minutes}))" \
            "$((10#${seconds}))"
        return
    fi
    if [[ "${slurm_time}" =~ ^([0-9]{1,3}):([0-9]{2}):([0-9]{2})$ ]]; then
        printf '%s\n' "${slurm_time}"
        return
    fi
    echo "Unsupported TIME_LIMIT format for PBS conversion: ${slurm_time}" >&2
    exit 1
}

ACTIVE_SCHEDULER=$(detect_scheduler)

if [[ "${ACTIVE_SCHEDULER}" == "slurm" ]]; then
    JOB_ID=$(
        sbatch --parsable \
            --partition="${PARTITION}" \
            --job-name="${JOB_NAME}" \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task="${CPUS_PER_TASK}" \
            --hint="${SLURM_HINT}" \
            --time="${TIME_LIMIT}" \
            --output="${SLURM_LOG_DIR}/%x-%j.out" \
            --error="${SLURM_LOG_DIR}/%x-%j.err" \
            --export=ALL,SCHEDULER_FAMILY=slurm,REPO_ROOT="${REPO_ROOT}",MEASUREMENT_CONFIG="${MEASUREMENT_CONFIG}",PYTHON_BIN="${PYTHON_BIN}",GCC_MODULE="${GCC_MODULE}",MONTE_CARLO_BUILD_DIR="${MONTE_CARLO_BUILD_DIR}",BATCH_SIZE="${BATCH_SIZE}",HIGHER_ORDER_BLOCK_SIZE="${HIGHER_ORDER_BLOCK_SIZE}" \
            "${SLURM_SCRIPT}"
    )
    JOB_ID_FILE="${SLURM_LOG_DIR}/${JOB_NAME}.latest_jobid"
    printf '%s\n' "${JOB_ID}" > "${JOB_ID_FILE}"

    echo "submitted_scheduler=slurm"
    echo "submitted_job_id=${JOB_ID}"
    echo "partition=${PARTITION}"
    echo "job_name=${JOB_NAME}"
    echo "time_limit=${TIME_LIMIT}"
    echo "cpus_per_task=${CPUS_PER_TASK}"
    echo "slurm_hint=${SLURM_HINT}"
    echo "monte_carlo_build_dir=${MONTE_CARLO_BUILD_DIR}"
    echo "batch_size=${BATCH_SIZE}"
    echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
    echo "resume_mode=enabled"
    echo "slurm_stdout=${SLURM_LOG_DIR}/${JOB_NAME}-${JOB_ID}.out"
    echo "slurm_stderr=${SLURM_LOG_DIR}/${JOB_NAME}-${JOB_ID}.err"
    echo "measurement_log=${RUN_LOG}"
    echo "progress_json=${PROGRESS_JSON}"
    echo "partial_rows_csv=${PARTIAL_ROWS}"
    echo "job_id_file=${JOB_ID_FILE}"
    exit 0
fi

if [[ ! -f "${PBS_SCRIPT}" ]]; then
    echo "Missing PBS wrapper: ${PBS_SCRIPT}" >&2
    exit 1
fi

SUBMIT_TAG=$(date -u +%Y%m%dT%H%M%SZ)
PBS_TIME_LIMIT=$(slurm_time_to_pbs "${TIME_LIMIT}")
PBS_STDOUT="${PBS_LOG_DIR}/${JOB_NAME}-${SUBMIT_TAG}.out"
PBS_STDERR="${PBS_LOG_DIR}/${JOB_NAME}-${SUBMIT_TAG}.err"
JOB_ID=$(
    qsub \
        -N "${JOB_NAME}" \
        -q "${QUEUE}" \
        -l "select=1:ncpus=${CPUS_PER_TASK}" \
        -l "walltime=${PBS_TIME_LIMIT}" \
        -o "${PBS_STDOUT}" \
        -e "${PBS_STDERR}" \
        -v "SCHEDULER_FAMILY=pbs,REPO_ROOT=${REPO_ROOT},MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG},PYTHON_BIN=${PYTHON_BIN},GCC_MODULE=${GCC_MODULE},MONTE_CARLO_BUILD_DIR=${MONTE_CARLO_BUILD_DIR},BATCH_SIZE=${BATCH_SIZE},HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE}" \
        "${PBS_SCRIPT}"
)
JOB_ID_FILE="${PBS_LOG_DIR}/${JOB_NAME}.latest_jobid"
printf '%s\n' "${JOB_ID}" > "${JOB_ID_FILE}"

echo "submitted_scheduler=pbs"
echo "submitted_job_id=${JOB_ID}"
echo "queue=${QUEUE}"
echo "job_name=${JOB_NAME}"
echo "time_limit=${PBS_TIME_LIMIT}"
echo "cpus_per_task=${CPUS_PER_TASK}"
echo "monte_carlo_build_dir=${MONTE_CARLO_BUILD_DIR}"
echo "batch_size=${BATCH_SIZE}"
echo "higher_order_block_size=${HIGHER_ORDER_BLOCK_SIZE}"
echo "resume_mode=enabled"
echo "pbs_stdout=${PBS_STDOUT}"
echo "pbs_stderr=${PBS_STDERR}"
echo "measurement_log=${RUN_LOG}"
echo "progress_json=${PROGRESS_JSON}"
echo "partial_rows_csv=${PARTIAL_ROWS}"
echo "job_id_file=${JOB_ID_FILE}"
