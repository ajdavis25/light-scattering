#!/usr/bin/env python3
"""Build the AERONET-constrained day-specific Marseille profile (2026-08-10).

Construction (user-confirmed "BL-weighted"): scale aerosol_extinction_550_m_inv
at grid points with altitude <= 700 m by a single factor s, chosen so the
trapezoid column AOD(550) equals the at-site AERONET value. Free troposphere
stays at the frozen climatology (0.087), which lies inside OHP's evening-to-
morning bracket 0.072-0.13. All other columns bit-identical to the frozen
profile. Aerosol optical properties (SSA/asymmetry/Angstrom) unchanged: same
fine-mode aerosol type (alpha ~ 1.2-1.5, no dust), just more of it in the BL.

Target provenance: Marseille_ATMO (43.3059 N, 5.3950 E, 65 m — at the
observation site; instrument #944) — the station whose direct-sun file was
ALREADY BUNDLED in the frozen case dir (aeronet_direct_sun.csv, L1.5); a fresh
V3 L2.0 pull is numerically identical (final calibration confirmed). Aug-15
afternoon-evening (14:00-17:06 UTC) plateau 0.1425-0.171, mean 0.155; last
point 0.1425 at 17:06; next-morning first 0.2167; linear bracket interpolation
to the 19:14 UTC observation ~0.155. TARGET = 0.155 (+-0.02 declared).
Regional context: Toulon 0.175 (upper bound), OHP FT-consistent.
"""
import csv

FROZEN = "/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/atmosphere_profile.csv"
OUT = "/work/vmo703/light-scattering/monte_carlo_cpp/data/atmosphere/marseille_profile_20220815_aeronet_bl.csv"
TARGET = 0.155   # Marseille_ATMO at-site evening value (see docstring)
BL_TOP = 700.0   # scale grid points at or below this altitude (m)

with open(FROZEN, newline="") as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = [row for row in reader]
icol = header.index("aerosol_extinction_550_m_inv")
zcol = header.index("altitude_m")

pts = [(float(r[zcol]), float(r[icol])) for r in rows]
assert pts == sorted(pts), "profile must be altitude-sorted"

def column(scale):
    vals = [e * scale if z <= BL_TOP else e for z, e in pts]
    return sum(0.5 * (vals[i] + vals[i + 1]) * (pts[i + 1][0] - pts[i][0])
               for i in range(len(pts) - 1))

# column(s) is linear in s: solve exactly from two evaluations
c0, c1 = column(1.0), column(2.0)
s = 1.0 + (TARGET - c0) / (c1 - c0)
achieved = column(s)
assert abs(achieved - TARGET) < 1e-12, achieved

bl_scaled = [i for i, (z, _) in enumerate(pts) if z <= BL_TOP]
print(f"frozen column = {c0:.6f}; scale s = {s:.6f} applied to grid points "
      f"{[pts[i][0] for i in bl_scaled]} m; achieved column = {achieved:.6f}")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(header)
    for r in rows:
        if float(r[zcol]) <= BL_TOP:
            r = list(r)
            r[icol] = repr(float(r[icol]) * s)
        w.writerow(r)
print(f"wrote {OUT}")
