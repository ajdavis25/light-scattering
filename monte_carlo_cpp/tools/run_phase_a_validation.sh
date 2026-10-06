#!/usr/bin/env bash
# Phase A validation rerun driver (see notebooks/PAPER_READINESS_PLAN_2026-07-26.md).
#
# Runs, sequentially and with a bounded thread count:
#   1. the default validation gate  (convergence + DISORT scalar + IPRT A1 vector
#      + Rozenberg + Koomen)                     -> default_clear_sky.cfg
#   2. the Zawada spherical-vector single-scatter smoke benchmark
#   3. the Zawada spherical-vector all-orders smoke benchmark
#
# ValidationRunner always writes monte_carlo_cpp/results/validation/validation_report.json
# (a fixed path), so this driver archives the report and the runner stdout under
# case-specific names after each run to prevent clobbering.
#
# Must be run from the repository root.
set -u

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

# RUNNER can be overridden (e.g. build_cluster_gcc8 for Slurm compute nodes).
RUNNER="${RUNNER:-monte_carlo_cpp/build_current/ValidationRunner}"
OUT_DIR=monte_carlo_cpp/results/validation
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
mkdir -p "$OUT_DIR"

echo "phase_a_driver_start utc=$STAMP omp_threads=$OMP_NUM_THREADS runner=$RUNNER"

run_case () {
    local label="$1"
    local cfg="$2"
    echo "=== run_case ${label} cfg=${cfg} start=$(date -u +%H:%M:%SZ) ==="
    nice -n 10 stdbuf -oL -eL "$RUNNER" "$cfg" > "$OUT_DIR/${label}_stdout.log" 2>&1
    local status=$?
    if [ -f "$OUT_DIR/validation_report.json" ]; then
        cp "$OUT_DIR/validation_report.json" "$OUT_DIR/validation_report__${label}.json"
    fi
    echo "=== run_case ${label} exit=${status} end=$(date -u +%H:%M:%SZ) ==="
    return $status
}

# Case list: label=cfg arguments, or the historical default list.
# Fastest-first so evidence lands early even if the job is cut short.
overall=0
if [ "$#" -gt 0 ]; then
    for pair in "$@"; do
        label="${pair%%=*}"
        cfg="${pair#*=}"
        run_case "$label" "$cfg" || overall=1
    done
else
    run_case zawada_single_smoke         monte_carlo_cpp/config/benchmark_zawada_spherical_vector_single_smoke.cfg || overall=1
    run_case zawada_multiple_smoke       monte_carlo_cpp/config/benchmark_zawada_spherical_vector_multiple_smoke.cfg || overall=1
    run_case default_clear_sky           monte_carlo_cpp/config/default_clear_sky.cfg || overall=1
fi

echo "phase_a_driver_done overall_exit=${overall} utc=$(date -u +%Y%m%dT%H%M%SZ)"
exit "$overall"
