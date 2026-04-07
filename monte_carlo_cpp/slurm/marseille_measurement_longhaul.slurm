#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=${REPO_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}
RUNNER_SCRIPT="${REPO_ROOT}/monte_carlo_cpp/tools/run_marseille_measurement_job.sh"

if [[ ! -f "${RUNNER_SCRIPT}" ]]; then
    echo "Missing runner script: ${RUNNER_SCRIPT}" >&2
    exit 1
fi

bash "${RUNNER_SCRIPT}"
