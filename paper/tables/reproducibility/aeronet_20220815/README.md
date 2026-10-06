# AERONET constraint on the 2022-08-15 Marseille aerosol loading (2026-08-10)

Question under test: does the frozen Marseille profile (column AOD(550 nm) =
0.1379) understate the real 2022-08-15 day by the ×2–×4 the DoLP loading scans
would require (0.28–0.55)? Observation: Marseille 43.28699 N, 5.40336 E,
2022-08-15 19:14:13 UTC (SZA 96.007).

## Verdict: ×2–×4 REFUTED; day-specific column ≈ 0.155 ± 0.02 (×1.12 frozen)

Three stations, all fine-mode (no dust; Ångström α ≈ 1.2–1.5):

| station | where | level | Aug-15 evening AOD(550) |
|---|---|---|---|
| **Marseille_ATMO** (#944) | AT the site (43.306 N, 5.395 E, 65 m) | 1.5 = 2.0 (bit-identical) | plateau 0.1425–0.171 (14:00–17:06 UTC), last point **0.1425** at 17:06; next morning 0.2167 |
| Toulon | coastal, 52 km ESE, 50 m | 1.5 | stable 0.175–0.184 through 17:46; next morning 0.1735 |
| OHP_OBSERVATOIRE | inland, ~75 km N, 680 m | 2.0 | declining to 0.072 by 17:49; next morning 0.124–0.155 |

- **At-site value at 19:14 UTC**: afternoon-evening mean 0.155; linear bracket
  interpolation (0.1425 at 17:06 → 0.2167 at 05:49+1d) also gives ≈ 0.155.
  Adopted target **0.155 ± 0.02**; Toulon 0.175 = regional upper bound.
- **Layering**: Toulon−OHP ≈ 0.10 below ~680 m vs frozen BL 0.0375 → the
  regional enhancement is boundary-layer; OHP's 0.072–0.13 bracket straddles
  the frozen free troposphere 0.087 → FT kept frozen in the v2 profile
  (BL-weighted scaling, s = 1.3796 below 700 m; `build_day_profile.py`).

## Provenance (important correction of the first-pass reading)

The frozen case dir ALREADY bundles the Marseille_ATMO direct-sun file
(`aeronet_direct_sun.csv`), and the v1 input pipeline was ALREADY day-informed:
`paper_case_provenance.json` records target_aod_550 = 0.1379 =
0.75 × AERONET(0.14393, anchored near the ~17:00 last direct-sun points)
+ 0.25 × Open-Meteo(0.12). The v1↔v2 difference (0.138 → 0.155, +12%) is an
evening-extrapolation choice on the SAME station — last-point anchor blended
with a model value, versus plateau-mean/bracket interpolation to 19:14 —
NOT a missing-data archaeology. A fresh V3 Level-2.0 pull of Marseille_ATMO is
numerically identical to the bundled L1.5 (final calibration confirmed).

## Implication for the DoLP gap

The ×2–×4 loading scans were a compensating knob, not a physical correction:
the licensed input-side day correction is ±10–15% in column. Measured on the
48-direction subset (v2 stack = chi-fix + King 0.0279 + 0.155 BL-weighted
profile): signed DoLP bias +0.1938 → +0.1678, of which the profile contributes
only ≈ −0.003 beyond King. **The residual ≈ +0.17 over-polarization is a
genuine open model/reduction discrepancy** (Marseille-vs-Koomen 4–6× bias
asymmetry = prime lead). Side effect worth keeping: the day profile alone
improves AoP median 23.4° → 18.0° (best recorded on the subset).

## Files

- `aeronet_marseille_atmo_L20.txt` — fresh V3 L2.0 pull, at-site station (PI
  from bundled file: Philippe Goloub / service processing; instrument #944)
- `aeronet_toulon_L15.txt` — Toulon (PI Jacques Piazzola)
- `aeronet_ohp_L20.txt` — OHP (PI Philippe Goloub)
- `parse_aeronet.py`, `aod550_summary.csv` — AOD(550) = AOD(500)·(1.1)^(−α),
  α = 440–675 Ångström exponent; per-day summaries
- `build_day_profile.py`, `marseille_profile_20220815_aeronet_bl.csv` — the v2
  day profile (column 0.155, BL-weighted, all non-aerosol columns bit-identical
  to frozen)
- Query URL form:
  `https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site=Marseille_ATMO&year=2022&month=8&day=15&year2=2022&month2=8&day2=16&AOD20=1&AVG=10&if_no_html=1`

Paper acknowledgment owed to the AERONET PIs of Marseille_ATMO, Toulon
(Jacques Piazzola), and OHP_OBSERVATOIRE (Philippe Goloub).
