#!/usr/bin/env python3
"""Score a measurement-case comparison CSV with the depol-scan summary metrics.

Reproduces the metric definitions of
paper/tables/reproducibility/marseille_depol_scan_summary_2026-08-10.csv:
  sgn_bias_med      median of signed_dolp_bias (model - reference DoLP)
  sgn_bias_zen_lt45 same, over rows with zenith < 45 deg
  aop_med_deg       median of aop_abs_error_deg
  shape_rmse        rmse(normalized_model - normalized_reference) (peak-normalized)
  model_dolp_med    median of model_dop
  abs_bias_med      median of dop_abs_error

Usage: score_measurement_comparison.py <comparison.csv> [label]
"""
import csv
import math
import sys


def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main():
    path = sys.argv[1]
    label = sys.argv[2] if len(sys.argv) > 2 else path
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    sgn = [float(r["signed_dolp_bias"]) for r in rows]
    sgn45 = [float(r["signed_dolp_bias"]) for r in rows if float(r["zenith_deg"]) < 45.0]
    aop = [float(r["aop_abs_error_deg"]) for r in rows]
    shape = [float(r["normalized_model"]) - float(r["normalized_reference"]) for r in rows]
    mdolp = [float(r["model_dop"]) for r in rows]
    absb = [float(r["dop_abs_error"]) for r in rows]
    rmse = math.sqrt(sum(d * d for d in shape) / len(shape))
    print(f"{label},n={len(rows)}")
    print(f"sgn_bias_med={median(sgn):+.4f} sgn_bias_zen_lt45={median(sgn45):+.4f} "
          f"aop_med_deg={median(aop):.2f} shape_rmse={rmse:.4f} "
          f"model_dolp_med={median(mdolp):.4f} abs_bias_med={median(absb):.4f}")


if __name__ == "__main__":
    main()
