#!/usr/bin/env python3
"""Parse AERONET v3 direct-sun files, interpolate AOD(550), summarize vs frozen profile.

AOD(550) = AOD(500) * (550/500)^(-alpha), alpha = 440-675 Angstrom exponent.
Frozen Marseille profile column AOD(550) = 0.1379; scan multiples x2=0.276 x3=0.414 x4=0.552.
Observation: 2022-08-15 19:14:13 UTC (evening twilight).
"""
import csv
import math
import sys

FROZEN = 0.1379

def parse(path):
    rows = []
    with open(path) as f:
        lines = f.read().splitlines()
    # find header line (starts with AERONET_Site,Date)
    hdr_i = next((i for i, l in enumerate(lines) if l.startswith("AERONET_Site,Date")), None)
    if hdr_i is None:
        return None, []
    hdr = lines[hdr_i].split(",")
    col = {name: i for i, name in enumerate(hdr)}
    i500 = col["AOD_500nm"]
    iae = col["440-675_Angstrom_Exponent"]
    isza = col["Solar_Zenith_Angle(Degrees)"]
    for l in lines[hdr_i + 1:]:
        p = l.split(",")
        if len(p) < isza + 1:
            continue
        try:
            a500 = float(p[i500]); ae = float(p[iae]); sza = float(p[isza])
        except ValueError:
            continue
        if a500 < 0 or ae < -900:
            continue
        d, t = p[1], p[2]  # dd:mm:yyyy, hh:mm:ss
        a550 = a500 * (550.0 / 500.0) ** (-ae)
        rows.append((d, t, a500, ae, a550, sza))
    return hdr, rows

def summarize(name, rows):
    print(f"\n===== {name}: {len(rows)} points =====")
    by_day = {}
    for r in rows:
        by_day.setdefault(r[0], []).append(r)
    for d in sorted(by_day, key=lambda s: s.split(":")[::-1]):
        pts = by_day[d]
        a550s = sorted(x[4] for x in pts)
        med = a550s[len(a550s) // 2]
        print(f"\n-- {d}: n={len(pts)} AOD550 min={a550s[0]:.4f} med={med:.4f} max={a550s[-1]:.4f} "
              f"(x{med/FROZEN:.2f} frozen)")
        for (dd, t, a500, ae, a550, sza) in pts:
            print(f"   {t}  AOD500={a500:.4f}  alpha={ae:.3f}  AOD550={a550:.4f}  SZA={sza:.1f}")

def main():
    for path, name in [("aeronet_toulon_L15.txt", "Toulon L1.5 (coastal)")]:
        try:
            hdr, rows = parse(path)
        except FileNotFoundError:
            print(f"{name}: file missing"); continue
        if hdr is None:
            print(f"{name}: no data header (empty return)"); continue
        summarize(name, rows)
    print(f"\nfrozen baseline AOD550 = {FROZEN}; x2 = {FROZEN*2:.3f}  x3 = {FROZEN*3:.3f}  x4 = {FROZEN*4:.3f}")

if __name__ == "__main__":
    main()
