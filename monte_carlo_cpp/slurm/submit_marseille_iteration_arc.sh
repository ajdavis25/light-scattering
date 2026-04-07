#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/../.." && pwd)
BASE_SUBMIT="${SCRIPT_DIR}/submit_marseille_measurement_arc.sh"

MODE=${MODE:-profile_subset}
PARTITION=${PARTITION:-compute1}
PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3.12}
GCC_MODULE=${GCC_MODULE:-gcc/14.2.0}

case "${MODE}" in
    profile_subset)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset.cfg"
        JOB_NAME_DEFAULT="marseille_profile_subset"
        CPUS_PER_TASK_DEFAULT=12
        BATCH_SIZE_DEFAULT=12
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=4
        TIME_LIMIT_DEFAULT="06:00:00"
        ;;
    strict_subset)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_strict_subset.cfg"
        JOB_NAME_DEFAULT="marseille_strict_subset"
        CPUS_PER_TASK_DEFAULT=48
        BATCH_SIZE_DEFAULT=48
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=4
        TIME_LIMIT_DEFAULT="12:00:00"
        ;;
    tiny_subset)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset.cfg"
        JOB_NAME_DEFAULT="marseille_tiny_subset"
        CPUS_PER_TASK_DEFAULT=16
        BATCH_SIZE_DEFAULT=16
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=2
        TIME_LIMIT_DEFAULT="02:00:00"
        ;;
    subset_quick)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_subset_quick.cfg"
        JOB_NAME_DEFAULT="marseille_subset_quick"
        CPUS_PER_TASK_DEFAULT=64
        BATCH_SIZE_DEFAULT=64
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=4
        TIME_LIMIT_DEFAULT="04:00:00"
        ;;
    smoke)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_smoke.cfg"
        JOB_NAME_DEFAULT="marseille_smoke"
        CPUS_PER_TASK_DEFAULT=80
        BATCH_SIZE_DEFAULT=80
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=8
        TIME_LIMIT_DEFAULT="08:00:00"
        ;;
    quick)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_quick.cfg"
        JOB_NAME_DEFAULT="marseille_quick"
        CPUS_PER_TASK_DEFAULT=80
        BATCH_SIZE_DEFAULT=80
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=4
        TIME_LIMIT_DEFAULT="08:00:00"
        ;;
    full_strict)
        MEASUREMENT_CONFIG_DEFAULT="${REPO_ROOT}/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg"
        JOB_NAME_DEFAULT="marseille_strict_single_node"
        CPUS_PER_TASK_DEFAULT=80
        BATCH_SIZE_DEFAULT=80
        HIGHER_ORDER_BLOCK_SIZE_DEFAULT=8
        TIME_LIMIT_DEFAULT="2-00:00:00"
        ;;
    *)
        echo "Unsupported MODE=${MODE}" >&2
        exit 1
        ;;
esac

export PARTITION
export PYTHON_BIN
export GCC_MODULE
export MEASUREMENT_CONFIG=${MEASUREMENT_CONFIG:-${MEASUREMENT_CONFIG_DEFAULT}}
export JOB_NAME=${JOB_NAME:-${JOB_NAME_DEFAULT}}
export CPUS_PER_TASK=${CPUS_PER_TASK:-${CPUS_PER_TASK_DEFAULT}}
export BATCH_SIZE=${BATCH_SIZE:-${BATCH_SIZE_DEFAULT}}
export HIGHER_ORDER_BLOCK_SIZE=${HIGHER_ORDER_BLOCK_SIZE:-${HIGHER_ORDER_BLOCK_SIZE_DEFAULT}}
export TIME_LIMIT=${TIME_LIMIT:-${TIME_LIMIT_DEFAULT}}

exec "${BASE_SUBMIT}"
