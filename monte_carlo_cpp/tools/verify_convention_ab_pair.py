"""Phase B5 A/B verification (PAPER_READINESS_PLAN_2026-07-26.md).

Compares the legacy vs conventionfix tiny-smoke comparison CSVs (same seed,
same binary) and checks the exact expected invariants of the incoming-light
convention transform (Q -> -Q):

  1. model_q(fix)  == -model_q(legacy)          (sign flip, exact)
  2. model_u(fix)  ==  model_u(legacy)          (unchanged, exact)
  3. model_dop     unchanged                    (invariant)
  4. model_intensity unchanged                  (invariant)
  5. model_aop(fix) == wrap(90 - model_aop(legacy))

Anything outside float tolerance fails loudly. Estimator note: the tiny smoke
config has no robust groups, so higher_order_estimator=mean is a no-op there
by design; the estimator branch is validated separately.

Usage (from repo root):
    python3 monte_carlo_cpp/tools/verify_convention_ab_pair.py
"""
import csv
import math
import os
import sys

REPORT_DIR = "monte_carlo_cpp/results/measurement_case_reports"
LEGACY = os.path.join(
    REPORT_DIR, "frozen_marseille_twilight_20220815_191413z_measurement_tiny_smoke_comparison.csv"
)
FIXED = os.path.join(
    REPORT_DIR,
    "frozen_marseille_twilight_20220815_191413z_measurement_tiny_smoke_conventionfix_comparison.csv",
)


def wrap_half_turn(v):
    while v <= -90.0:
        v += 180.0
    while v > 90.0:
        v -= 180.0
    return v


def load(path):
    with open(path, newline="") as handle:
        return {int(float(r["index"])): r for r in csv.DictReader(handle)}


def main():
    legacy = load(LEGACY)
    fixed = load(FIXED)
    if set(legacy) != set(fixed):
        print("FAIL: index sets differ: %d vs %d" % (len(legacy), len(fixed)))
        return 1

    tol = 1e-9
    failures = 0
    for idx in sorted(legacy):
        a, b = legacy[idx], fixed[idx]
        lq, lu = float(a["model_q"]), float(a["model_u"])
        fq, fu = float(b["model_q"]), float(b["model_u"])
        li, fi = float(a["model_intensity"]), float(b["model_intensity"])
        ld, fd = float(a["model_dop"]), float(b["model_dop"])
        la, fa = float(a["model_aop_deg"]), float(b["model_aop_deg"])

        scale = max(abs(lq), abs(lu), abs(li) * 1e-3, 1e-300)
        checks = [
            ("q_sign_flip", abs(fq + lq) / scale),
            ("u_unchanged", abs(fu - lu) / scale),
            ("intensity_unchanged", abs(fi - li) / max(abs(li), 1e-300)),
            ("dolp_unchanged", abs(fd - ld) / max(abs(ld), 1e-12)),
        ]
        aop_expected = wrap_half_turn(90.0 - la)
        aop_diff = abs(wrap_half_turn(fa - aop_expected))
        aop_diff = min(aop_diff, abs(180.0 - aop_diff))
        checks.append(("aop_90_minus", aop_diff / 90.0))

        for name, rel in checks:
            if rel > tol:
                print("FAIL idx=%d %s rel=%.3e" % (idx, name, rel))
                failures += 1

    if failures == 0:
        print("PASS: all %d directions satisfy the incoming-light transform invariants exactly."
              % len(legacy))
        return 0
    print("FAILURES: %d" % failures)
    return 1


if __name__ == "__main__":
    sys.exit(main())
