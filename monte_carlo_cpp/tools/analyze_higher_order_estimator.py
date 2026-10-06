"""Phase C1 (PAPER_READINESS_PLAN_2026-07-26.md): median-of-means vs plain mean.

The frozen full-field Marseille run (robust_r2) recorded, per direction, BOTH
estimators' ingredients in its checkpoint files:

  - higher_mean_{I,Q,U,V}: the plain running mean over all completed
    higher-order samples, and
  - higher_robust_group_{k}_{samples,sum_*}: the 16 robust-group sums used by
    RobustHigherOrderGroups::medianOfMeans, which is what the solver actually
    reported (MonteCarloDriver.cpp:3333).

So the estimator-bias question can be answered exactly, for the exact frozen
field, with no solver reruns. This script parses all checkpoints, reconstructs
both per-direction estimates, verifies the median-of-means reconstruction
against the raw model values in the frozen comparison CSV, and quantifies the
difference in higher-order I and in total-field I / DoLP / AoP.

Usage (from repo root):
    python3 monte_carlo_cpp/tools/analyze_higher_order_estimator.py
"""
import csv
import glob
import math
import os
import re
import statistics

REPORT_DIR = "monte_carlo_cpp/results/measurement_case_reports"
CHECKPOINT_GLOB = os.path.join(
    REPORT_DIR, "_batched_work", "*robust_r2*checkpoint.txt"
)
COMPARISON_CSV = os.path.join(
    REPORT_DIR,
    "frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2_comparison.csv",
)
OUTPUT_CSV = os.path.join(REPORT_DIR, "higher_order_estimator_comparison.csv")


def median_component(values):
    """Exact replication of medianComponent in MonteCarloDriver.cpp:100-113."""
    if not values:
        return 0.0
    ordered = sorted(values)
    middle = len(values) // 2
    median = ordered[middle]
    if len(values) % 2 == 0:
        median = 0.5 * (ordered[middle - 1] + median)
    return median


def parse_checkpoint(path):
    data = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            data[key] = value

    out = {
        "direction_index": int(data["direction_index"]),
        "zenith_deg": float(data["zenith_deg"]),
        "azimuth_deg": float(data["azimuth_deg"]),
        "samples": int(data.get("higher_completed_samples", "0")),
    }
    for comp in "IQUV":
        out[f"first_{comp}"] = float(data.get(f"first_order_{comp}", "0"))
        out[f"second_{comp}"] = float(data.get(f"second_total_{comp}", "0"))
        out[f"mean_{comp}"] = float(data.get(f"higher_mean_{comp}", "0"))

    group_count = int(data.get("higher_robust_group_count", "0"))
    group_means = {c: [] for c in "IQUV"}
    used_groups = 0
    for k in range(group_count):
        n = int(data.get(f"higher_robust_group_{k}_samples", "0"))
        if n <= 0:
            continue
        used_groups += 1
        for comp in "IQUV":
            total = float(data.get(f"higher_robust_group_{k}_sum_{comp}", "0"))
            group_means[comp].append(total / n)

    if used_groups >= 2:
        for comp in "IQUV":
            out[f"mom_{comp}"] = median_component(group_means[comp])
    else:  # medianOfMeans falls back to the plain mean
        for comp in "IQUV":
            out[f"mom_{comp}"] = out[f"mean_{comp}"]
    out["used_groups"] = used_groups
    out["group_mean_I_spread"] = (
        max(group_means["I"]) / max(min(group_means["I"]), 1e-300)
        if group_means["I"] and min(group_means["I"]) > 0
        else float("nan")
    )
    return out


def dolp(i, q, u):
    return math.hypot(q, u) / i if i > 0 else 0.0


def aop_deg(q, u):
    return 0.5 * math.degrees(math.atan2(u, q))


def wrap_half_turn(value):
    while value <= -90.0:
        value += 180.0
    while value > 90.0:
        value -= 180.0
    return value


def pct(values, q):
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    pos = min(len(ordered) - 1, max(0, int(round(q * (len(ordered) - 1)))))
    return ordered[pos]


def main():
    paths = sorted(glob.glob(CHECKPOINT_GLOB))
    print(f"checkpoints_found={len(paths)}")
    rows = [parse_checkpoint(p) for p in paths]
    rows.sort(key=lambda r: r["direction_index"])

    # Cross-check: does first + second + median-of-means reproduce the raw
    # model values actually scored in the frozen comparison CSV?
    raw_by_index = {}
    with open(COMPARISON_CSV, "r", encoding="utf-8", newline="") as handle:
        for rec in csv.DictReader(handle):
            raw_by_index[int(float(rec["index"]))] = rec

    reconstruction_mismatches = 0
    worst_recon_rel = 0.0
    for r in rows:
        rec = raw_by_index.get(r["direction_index"])
        if rec is None:
            continue
        recon_i = r["first_I"] + r["second_I"] + r["mom_I"]
        ref_raw_i = float(rec["raw_model_intensity"])
        if ref_raw_i > 0:
            rel = abs(recon_i - ref_raw_i) / ref_raw_i
            worst_recon_rel = max(worst_recon_rel, rel)
            if rel > 1e-9:
                reconstruction_mismatches += 1
    print(
        f"reconstruction_check: mismatches(>1e-9 rel)={reconstruction_mismatches}"
        f"/{len(rows)} worst_rel={worst_recon_rel:.3e}"
    )

    # Per-direction estimator deltas.
    higher_rel_diff = []       # (mom - mean)/mean on higher-order I
    total_rel_diff = []        # same on total I
    dolp_delta = []            # DoLP(mom-total) - DoLP(mean-total)
    aop_delta = []             # wrapped AoP difference, deg
    higher_frac_mom = []
    total_i_mom_sum = 0.0
    total_i_mean_sum = 0.0

    out_rows = []
    for r in rows:
        tm_i = r["first_I"] + r["second_I"] + r["mean_I"]
        tm_q = r["first_Q"] + r["second_Q"] + r["mean_Q"]
        tm_u = r["first_U"] + r["second_U"] + r["mean_U"]
        to_i = r["first_I"] + r["second_I"] + r["mom_I"]
        to_q = r["first_Q"] + r["second_Q"] + r["mom_Q"]
        to_u = r["first_U"] + r["second_U"] + r["mom_U"]

        total_i_mom_sum += to_i
        total_i_mean_sum += tm_i

        h_rel = (r["mom_I"] - r["mean_I"]) / r["mean_I"] if r["mean_I"] > 0 else 0.0
        t_rel = (to_i - tm_i) / tm_i if tm_i > 0 else 0.0
        d_delta = dolp(to_i, to_q, to_u) - dolp(tm_i, tm_q, tm_u)
        a_delta = wrap_half_turn(aop_deg(to_q, to_u) - aop_deg(tm_q, tm_u))

        higher_rel_diff.append(h_rel)
        total_rel_diff.append(t_rel)
        dolp_delta.append(d_delta)
        aop_delta.append(a_delta)
        higher_frac_mom.append(r["mom_I"] / to_i if to_i > 0 else 0.0)

        out_rows.append({
            "direction_index": r["direction_index"],
            "zenith_deg": r["zenith_deg"],
            "azimuth_deg": r["azimuth_deg"],
            "samples": r["samples"],
            "used_groups": r["used_groups"],
            "higher_I_mean": f"{r['mean_I']:.9e}",
            "higher_I_mom": f"{r['mom_I']:.9e}",
            "higher_I_rel_diff": f"{h_rel:.6e}",
            "group_mean_I_spread": f"{r['group_mean_I_spread']:.3e}",
            "total_I_rel_diff": f"{t_rel:.6e}",
            "dolp_delta": f"{d_delta:.6e}",
            "aop_delta_deg": f"{a_delta:.6e}",
            "higher_frac_mom": f"{higher_frac_mom[-1]:.6e}",
        })

    with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"wrote {OUTPUT_CSV} rows={len(out_rows)}")

    def summarize(name, values, absolute=True):
        source = [abs(v) for v in values] if absolute else list(values)
        print(
            f"{name}: median={statistics.median(source):.4e} "
            f"p95={pct(source, 0.95):.4e} max={max(source):.4e} "
            f"mean_signed={statistics.mean(values):+.4e}"
        )

    print("\n--- higher-order component: (median-of-means - mean)/mean on I ---")
    summarize("higher_I_rel_diff", higher_rel_diff)
    negative = sum(1 for v in higher_rel_diff if v < 0)
    print(f"directions where median-of-means < plain mean: {negative}/{len(rows)}")

    print("\n--- impact on the TOTAL per-direction field ---")
    summarize("total_I_rel_diff", total_rel_diff)
    summarize("dolp_delta", dolp_delta)
    summarize("aop_delta_deg", aop_delta)
    print(
        "field-wide summed intensity: mom/mean = "
        f"{total_i_mom_sum / total_i_mean_sum:.6f}"
    )
    print(
        f"median higher-order fraction of total I (mom): "
        f"{statistics.median(higher_frac_mom):.4f}"
    )

    # Context: compare against the raw Marseille misfit scale so the estimator
    # question can be judged as material or immaterial.
    print(
        "\ncontext: raw Marseille misfit is ~0.19 median |dDoLP| and ~44 deg "
        "mean |dAoP| (audit M-1); compare magnitudes above against that."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
