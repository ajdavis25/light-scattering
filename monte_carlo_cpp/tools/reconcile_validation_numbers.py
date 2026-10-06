"""Phase A5 (PAPER_READINESS_PLAN_2026-07-26.md): reconcile regenerated
validation metrics against the values historically quoted in STATUS_README.md,
MODEL_ASSUMPTIONS.md, and monte_carlo_cpp/README.md.

Reads the per-case archived reports produced by tools/run_phase_a_validation.sh
(validation_report__<label>.json) and prints fresh-vs-claimed side by side.

Each claim family lists the report labels that may carry it, tried in order:
the four benchmark/measurement families can come from either the split
single-case wrappers (validation_split_*.cfg, cheap, land first) or the full
default gate; the convergence family may ONLY come from the full default gate
(the split wrappers shrink the top-level sky to 1x1 bins at one wavelength, so
their convergence_* metrics are meaningless by design).

Exact bit-for-bit agreement is NOT expected: the historical numbers came from
the Windows MinGW build with different thread counts (OpenMP reduction order
changes float accumulation). Same-threshold pass/fail and same order of
magnitude of each metric is the reconciliation standard.

Usage (from repo root):
    python3 monte_carlo_cpp/tools/reconcile_validation_numbers.py
"""
import json
import os
import sys

VAL_DIR = "monte_carlo_cpp/results/validation"

# family -> (candidate report labels in preference order,
#            {metric name: (historically claimed value, where quoted)})
CLAIM_FAMILIES = [
    # default_clear_sky_capped = the 2026-08-03 declared-deviation gate
    # (validation_default_capped.cfg, higher_order_recursive_branch_cap=1):
    # the recorded uncapped config is non-terminating at its SZA-97 sun, so
    # the capped report is the only obtainable fresh convergence family.
    ("convergence", ["default_clear_sky", "default_clear_sky_capped"], {
        "convergence_peak_intensity_rel": (0.0382187, "STATUS_README.md"),
        "convergence_peak_dolp_abs": (0.00138362, "STATUS_README.md"),
        "convergence_flux_rel": (0.00211343, "STATUS_README.md"),
    }),
    ("disort_scalar", ["split_disort", "default_clear_sky"], {
        "benchmark_benchmark_disort_scalar_median_intensity_rel": (0.0188119, "STATUS_README.md"),
        "benchmark_benchmark_disort_scalar_p95_intensity_rel": (0.0250578, "STATUS_README.md"),
    }),
    ("iprt_a1_vector", ["split_iprt", "default_clear_sky"], {
        "benchmark_benchmark_iprt_a1_vector_median_intensity_rel": (0.00561889, "STATUS_README.md"),
        "benchmark_benchmark_iprt_a1_vector_p95_intensity_rel": (0.0176791, "STATUS_README.md"),
        "benchmark_benchmark_iprt_a1_vector_median_dolp_abs": (0.0055043, "STATUS_README.md"),
        "benchmark_benchmark_iprt_a1_vector_p95_dolp_abs": (0.0286531, "STATUS_README.md"),
    }),
    ("rozenberg", ["split_rozenberg", "default_clear_sky"], {
        "measurement_measurement_rozenberg_hminus6_normalized_rmse": (0.0849321, "STATUS_README.md"),
    }),
    ("koomen", ["split_koomen", "default_clear_sky"], {
        "measurement_measurement_koomen_meridian_hminus6_polarization_normalized_rmse": (0.0470789, "STATUS_README.md"),
        "measurement_measurement_koomen_meridian_hminus6_polarization_p95_dolp_abs": (0.0366844, "STATUS_README.md"),
    }),
    # split_zawada_base_legacy = the base-tier zawada cases (the tier the
    # claims came from); the *_smoke tier is stale and mismatches the bundled
    # reference, kept only as a fallback label for historical reports.
    ("zawada_single", ["split_zawada_base_legacy", "zawada_single_smoke"], {
        "benchmark_benchmark_zawada_spherical_vector_single_median_intensity_rel": (0.000459695, "STATUS_README/README"),
        "benchmark_benchmark_zawada_spherical_vector_single_p95_intensity_rel": (0.00140769, "STATUS_README/README"),
        "benchmark_benchmark_zawada_spherical_vector_single_median_dolp_abs": (3.37895e-05, "STATUS_README/README"),
        "benchmark_benchmark_zawada_spherical_vector_single_p95_dolp_abs": (0.000292954, "STATUS_README/README"),
    }),
    ("zawada_multiple", ["split_zawada_base_legacy", "zawada_multiple_smoke"], {
        "benchmark_benchmark_zawada_spherical_vector_multiple_median_intensity_rel": (0.00457567, "STATUS_README/README"),
        "benchmark_benchmark_zawada_spherical_vector_multiple_p95_intensity_rel": (0.012131, "STATUS_README/README"),
        "benchmark_benchmark_zawada_spherical_vector_multiple_median_dolp_abs": (0.0199496, "STATUS_README/README"),
        "benchmark_benchmark_zawada_spherical_vector_multiple_p95_dolp_abs": (0.0353911, "STATUS_README/README"),
    }),
]


def load_report(label):
    path = os.path.join(VAL_DIR, "validation_report__%s.json" % label)
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def main():
    families_missing = 0
    for family, labels, claims in CLAIM_FAMILIES:
        report = None
        used_label = None
        for label in labels:
            report = load_report(label)
            if report is not None:
                used_label = label
                break
        print("\n=== %s ===" % family)
        if report is None:
            print("  (no report yet; looked for labels: %s)" % ", ".join(labels))
            families_missing += 1
            continue
        print("  from validation_report__%s.json (overall_pass=%s)"
              % (used_label, report.get("overall_pass")))
        fresh = {m["name"]: m for m in report.get("metrics", [])}
        for name, (claimed, source) in sorted(claims.items()):
            if name in fresh:
                m = fresh[name]
                ratio = m["value"] / claimed if claimed else float("nan")
                print("  %-70s fresh=%.6g claimed=%.6g ratio=%.3f pass=%s"
                      % (name, m["value"], claimed, ratio, m["pass"]))
            else:
                print("  %-70s MISSING from fresh report (claimed %.6g in %s)"
                      % (name, claimed, source))
    print("\nfamilies missing: %d/%d" % (families_missing, len(CLAIM_FAMILIES)))
    return 1 if families_missing == len(CLAIM_FAMILIES) else 0


if __name__ == "__main__":
    sys.exit(main())
