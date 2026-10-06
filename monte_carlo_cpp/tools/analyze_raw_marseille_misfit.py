"""Phase B1 (PAPER_READINESS_PLAN_2026-07-26.md): decompose the raw Marseille misfit.

Inputs: the frozen calibration CSV, which stores per-direction raw model and
reference values (intensity, DoLP, AoP) for all 683 scored sky bins, plus the
frozen comparison CSV for cross-checks.

Questions answered here (python 3.6 compatible, numpy only):
  1. Absolute intensity scale: the measurement gate normalizes each field by
     its own peak (measurement_case_main.cpp:1635-1638), so absolute scale is
     invisible to the gate. Quantify the DN-vs-radiance scale factor anyway,
     and how much of the per-point gain spread survives after removing one
     global constant.
  2. Linearity: log-log slope of reference vs model intensity. Slope 1 means a
     pure scale factor; slope != 1 means a nonlinear response (exposure/
     gamma/saturation) or a systematic brightness-dependent model error.
  3. Shape: peak-normalized intensity residuals, masked like the gate
     (>= 5% of reference peak), and their correlation with viewing geometry.
  4. Where the DoLP/AoP failures live: vs zenith, vs scattering angle from the
     sun, and vs sky region.

Usage (from repo root):
    python3 monte_carlo_cpp/tools/analyze_raw_marseille_misfit.py
"""
import csv
import math
import os

import numpy as np

CASE_DIR = "monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z"
CAL_CSV = os.path.join(CASE_DIR, "measurement_model_quality_calibration.csv")
SUN_ZENITH_DEG = 96.007067  # from the frozen measurement cfg
OUT_CSV = "monte_carlo_cpp/results/measurement_case_reports/raw_misfit_decomposition.csv"


def wrap_half_turn(v):
    while v <= -90.0:
        v += 180.0
    while v > 90.0:
        v -= 180.0
    return v


def corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def main():
    rows = list(csv.DictReader(open(CAL_CSV)))
    n = len(rows)
    zen = np.array([float(r["zenith_deg"]) for r in rows])
    relaz = np.array([float(r["relative_azimuth_deg"]) for r in rows])
    ref_i = np.array([float(r["source_reference_intensity"]) for r in rows])
    mod_i = np.array([float(r["source_model_intensity"]) for r in rows])
    ref_d = np.array([float(r["source_reference_dop"]) for r in rows])
    mod_d = np.array([float(r["source_model_dop"]) for r in rows])
    ref_a = np.array([float(r["source_reference_aop_deg"]) for r in rows])
    mod_a = np.array([float(r["source_model_aop_deg"]) for r in rows])
    print("rows=%d" % n)

    # Scattering angle from the sun for each viewing direction.
    sz = math.radians(SUN_ZENITH_DEG)
    vz = np.radians(zen)
    dphi = np.radians(relaz)
    cos_theta = np.cos(vz) * math.cos(sz) + np.sin(vz) * math.sin(sz) * np.cos(dphi)
    scat_deg = np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))

    # ---------- 1. absolute scale + 2. linearity ----------
    log_ref = np.log10(ref_i)
    log_mod = np.log10(mod_i)
    slope, intercept = np.polyfit(log_mod, log_ref, 1)
    pred = slope * log_mod + intercept
    resid_loglog = log_ref - pred
    r2_slope = 1.0 - resid_loglog.var() / log_ref.var()

    # Pure-scale model (slope forced to 1):
    log_gain = log_ref - log_mod
    const_gain = log_gain.mean()          # best single constant, log space
    resid_const = log_gain - const_gain
    r2_const = 1.0 - resid_const.var() / log_ref.var()

    print("\n--- absolute scale & linearity (log10 space) ---")
    print("best single constant gain: 10^%.4f = %.3e" % (const_gain, 10 ** const_gain))
    print("residual spread after ONE constant: stdev=%.4f dex (x%.2f at 1-sigma)"
          % (resid_const.std(), 10 ** resid_const.std()))
    print("free log-log fit: slope=%.4f (1.0 = pure scale), intercept=%.3f" % (slope, intercept))
    print("variance explained: const-only R2=%.4f, free-slope R2=%.4f" % (r2_const, r2_slope))
    print("residual (const model) correlations: zenith=%.3f relaz=%.3f scat_angle=%.3f"
          % (corr(resid_const, zen), corr(resid_const, relaz), corr(resid_const, scat_deg)))
    print("residual (free-slope model) correlations: zenith=%.3f scat_angle=%.3f"
          % (corr(resid_loglog, zen), corr(resid_loglog, scat_deg)))

    # Reference DN distribution vs 12-bit range (saturation/compression check).
    print("reference intensity (DN?): min=%.1f p50=%.1f p95=%.1f max=%.1f  (12-bit max=4095)"
          % (ref_i.min(), np.percentile(ref_i, 50), np.percentile(ref_i, 95), ref_i.max()))

    # ---------- 3. gate-equivalent shape error ----------
    nref = ref_i / ref_i.max()
    nmod = mod_i / mod_i.max()
    mask = nref >= 0.05  # measurement_mask_fraction_of_peak
    shape_resid = nmod - nref
    rmse_all = float(np.sqrt((shape_resid ** 2).mean()))
    rmse_mask = float(np.sqrt((shape_resid[mask] ** 2).mean()))
    print("\n--- peak-normalized intensity shape error (the thing the gate actually scores) ---")
    print("masked points (>=5%% of ref peak): %d/%d" % (int(mask.sum()), n))
    print("normalized RMSE: all=%.4f masked=%.4f  (raw interactive runs reported 0.12-0.16)"
          % (rmse_all, rmse_mask))
    print("shape residual correlations (masked): zenith=%.3f scat_angle=%.3f"
          % (corr(shape_resid[mask], zen[mask]), corr(shape_resid[mask], scat_deg[mask])))
    # Where is the model too bright / too dim?
    bright = shape_resid[mask] > 0
    print("masked bins where model (normalized) > reference: %d/%d" % (int(bright.sum()), int(mask.sum())))

    # ---------- 4. DoLP / AoP failure structure ----------
    d_err = mod_d - ref_d
    a_err = np.array([wrap_half_turn(m - r) for m, r in zip(mod_a, ref_a)])
    print("\n--- raw DoLP misfit structure ---")
    print("signed DoLP error: mean=%+.4f median=%+.4f  (positive = model too polarized)"
          % (d_err.mean(), float(np.median(d_err))))
    print("|DoLP error|: median=%.4f p95=%.4f" % (float(np.median(np.abs(d_err))), np.percentile(np.abs(d_err), 95)))
    print("model DoLP range: %.3f-%.3f   reference DoLP range: %.3f-%.3f"
          % (mod_d.min(), mod_d.max(), ref_d.min(), ref_d.max()))
    print("signed error correlations: zenith=%.3f scat_angle=%.3f"
          % (corr(d_err, zen), corr(d_err, scat_deg)))
    # DoLP vs the Rayleigh single-scatter expectation shape:
    rayleigh_dolp = (1 - cos_theta ** 2) / (1 + cos_theta ** 2)
    print("corr(model DoLP, Rayleigh ss shape)=%.3f   corr(ref DoLP, Rayleigh ss shape)=%.3f"
          % (corr(mod_d, rayleigh_dolp), corr(ref_d, rayleigh_dolp)))

    print("\n--- raw AoP misfit structure ---")
    print("|AoP error|: median=%.2f deg p95=%.2f deg" % (float(np.median(np.abs(a_err))), np.percentile(np.abs(a_err), 95)))
    for lo, hi in [(0, 30), (30, 60), (60, 90)]:
        sel = (zen >= lo) & (zen < hi)
        if sel.sum():
            print("  zenith %2d-%2d deg (n=%3d): median |AoP err|=%.2f deg, median |DoLP err|=%.3f"
                  % (lo, hi, int(sel.sum()), float(np.median(np.abs(a_err[sel]))), float(np.median(np.abs(d_err[sel])))))
    # Is the AoP error consistent with a global basis rotation?
    print("signed AoP error: circular-mean=%.2f deg, stdev=%.2f deg"
          % (float(np.degrees(0.5 * np.angle(np.exp(2j * np.radians(a_err)).mean()))), a_err.std()))
    # A pure reference-frame error would give a tight distribution around a constant.

    # Save per-row decomposition for downstream plotting.
    with open(OUT_CSV, "w", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(["index", "zenith_deg", "relative_azimuth_deg", "scattering_angle_deg",
                    "log10_gain", "resid_const_dex", "norm_shape_resid", "masked",
                    "dolp_err_signed", "aop_err_signed_deg"])
        for k in range(n):
            w.writerow([rows[k]["index"], "%.4f" % zen[k], "%.4f" % relaz[k], "%.4f" % scat_deg[k],
                        "%.6f" % log_gain[k], "%.6f" % resid_const[k], "%.6f" % shape_resid[k],
                        int(mask[k]), "%.6f" % d_err[k], "%.4f" % a_err[k]])
    print("\nwrote %s" % OUT_CSV)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
