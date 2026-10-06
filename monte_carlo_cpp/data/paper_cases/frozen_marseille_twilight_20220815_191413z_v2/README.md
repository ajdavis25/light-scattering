# Frozen Marseille twilight tier V2 (2026-08-11, task #23)

Physics-only re-freeze of the 2022-08-15 19:14:13 UTC Marseille twilight case
(SZA 96.007). The v1 tier (`../frozen_marseille_twilight_20220815_191413z/`)
is UNTOUCHED; every input not listed below is inherited from it (measurement
reference, instrument response, optics tables, aerosol phase matrix, surface).

## The v2 stack — every non-default traces to physics or measurement

| knob | value | why |
|---|---|---|
| `event_frame_chi_sign_fix` | true | R(−χin)/R(+χout) handedness inconsistency, proven by dipole double-scatter ground truth (2026-08-01); vindicates IPRT claims to 4+ digits |
| `rayleigh_depolarization_factor` | 0.0279 | molecular (King-factor) anisotropy, Bates 1984; DoLP(90°)=0.9457 |
| `atmosphere_profile.csv` | column AOD(550)=0.155 | at-site AERONET (Marseille_ATMO) evening plateau/bracket at observation time; BL-weighted (s=1.3796 below 700 m, FT frozen — inside OHP bracket, Toulon−OHP layering); md5 bb7c477b74983fdcdf5eb837249b3d12 |
| `aerosol_depolarization_f22_ratio` | 1.0 (EXCLUDED) | no dust that day (α≈1.2–1.5 all three stations); f22=0.7 remains a documented sensitivity variant (AoP −3°, DoLP −0.008) |
| `higher_order_recursive_branch_cap` | 1 | declared deviation carried from v1 (uncapped recursion non-terminating, ~2^61) |

v1→v2 profile note: v1's 0.1379 was ALREADY day-informed
(0.75×AERONET(0.14393)+0.25×OpenMeteo(0.12), anchored on the same station's
last ~17:00 direct-sun points; see v1 `paper_case_provenance.json`). v2's
0.155 is a refinement of the same at-site data (plateau mean + evening-bracket
interpolation to 19:14 UTC), NOT a missing-data correction. External pulls +
verdict: `paper/tables/reproducibility/aeronet_20220815/`.

## Frozen results (jobs 788181 subset+full, 788184 gate; pinned binary build_v2freeze)

Comparison-CSV metric family (`tools/score_measurement_comparison.py`,
calibrated exact against the 2026-08-10 scan table):

| case | n | sgn_bias_med | aop_med_deg | shape_rmse | abs_bias_med |
|---|---|---|---|---|---|
| subset chifix baseline | 48 | +0.1938 | 23.37 | 0.1731 | 0.2288 |
| **subset v2** | 48 | **+0.1678** | **18.01** | 0.2599 | 0.1806 |
| full chifix baseline | 683 | +0.2068 | 17.01 | 0.0827 | 0.2173 |
| **full v2** | 683 | **+0.1736** | **16.96** | 0.0940 | 0.1981 |

Full-field zenith-bucket signed-bias medians (0-30/30-60/60-90°):
chifix +0.350/+0.143/+0.093 → v2 +0.304/+0.106/+0.076 (MS-dilution
fingerprint persists, all buckets improved).

**Declared open discrepancy:** the residual ≈ +0.17 median over-polarization
is NOT input-side (the licensed loading correction is ±10–15% and contributes
only ≈ −0.003 beyond King). Prime lead: Marseille-vs-Koomen 4–6× bias
asymmetry (reduction-side suspect) — Koomen fully passes under this identical
physics stack.

## V2 validation gate (validation_default_capped_v2.cfg, report archived)

22/23 metrics PASS (`validation_report__default_clear_sky_capped_v2.json`):
- IPRT A1 chifix: median_dolp_abs 0.0055043 (vindicated historical value)
- Rozenberg (true h⁻⁶ profile + King): normalized_rmse **0.0851 ≈ claim 0.0849**
- **Koomen (same h⁻⁶ profile + King): FULL PASS, first time** — RMSE 0.040
  (was ~0.10), median_dolp_abs 0.0160 < 0.02 (was 0.0515). Same input
  archaeology as Rozenberg, applied coherently to the same scenario family.
- Convergence: hardened metric (`convergence_low_order_metric=true` — peaks
  from deterministic first+second order; flux stays total-field), photon
  budget 1024.
- **Declared deviation (the 1 fail):** `convergence_flux_rel` = 0.0247 vs
  thr 0.02 — the characterized heavy-tailed higher-order estimator; PROVEN
  budget-independent (0.0240 at 256 photons, run 788183; 0.0247 at 1024, run
  788184). The legacy chain's 0.006 was an artifact of higher-order
  suppression (median_of_means + chi misrotation).

## Reproduction

```
export MONTE_CARLO_BUILD_DIR=.../monte_carlo_cpp/build_v2freeze   # PINNED — never rebuild
python3.12 monte_carlo_cpp/tools/run_measurement_case_batched.py \
    monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_full_v2.cfg \
    --batch-size 32 --higher-order-block-size 1 --resume
RUNNER=.../build_v2freeze/ValidationRunner bash monte_carlo_cpp/tools/run_phase_a_validation.sh \
    default_clear_sky_capped_v2=monte_carlo_cpp/config/validation_default_capped_v2.cfg
```

Binary md5: MeasurementCaseRunner 7391e37ec6d411c5d5b9f09af319874e,
ValidationRunner c2d3afca7dfd7ff7688494dcf13e50d8 (= build_fix3 source state,
2026-08-11; all-default outputs bit-identical to build_fix2, ctest 5/5).
RNG is case_id-seeded → reruns are bit-identical. Config copies in `configs/`.
Operational rule learned at cost (job 788177, ETXTBSY): NEVER rebuild a build
dir that a queued/running job executes from — dev work goes in build_fix3.
