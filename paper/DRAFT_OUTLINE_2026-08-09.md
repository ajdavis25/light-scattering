# Paper outline + skeleton (Phase E draft, 2026-08-09)

Status: DRAFT for the papers's direction decision. This file does not modify the
frozen April package (README/CLAIMS/MANIFEST/provenance stay authoritative for
the frozen artifacts); it lays out the candidate paper built on everything
learned July 26 - August 9. Canonical evidence log:
`notebooks/PLAN_IMPLEMENTATION_PROGRESS_2026-07-26.md`;
archived numbers: `paper/tables/reproducibility/`.

## Two candidate framings

**A. Calibrated-pipeline validation paper (the frozen April framing).**
The defensible-but-modest claim already packaged here: reproducible polarized
twilight RT pipeline + benchmark checks + calibrated full-field Marseille
artifact. Weakness: the calibrated pass validates the data path, not the
physics; referees could ask what the raw model does, and the honest answer
(raw twilight DoLP biased +0.19) undercuts the paper unless it is itself the
subject.

**B. Model + validation + regression-archaeology paper (recommended).**
A backward polarized Monte Carlo sky model with a validation suite treated as
a first-class scientific object: reproduced families, an exactly-diagnosed
polarization regression (found by independent dipole ground truth, vindicating
historical claims to 4+ significant figures), declared deviations where the
recorded configurations are non-terminating (with a 9-day/21-second
control pair), and an honestly quantified open twilight depolarization gap.
This converts every awkward reconciliation fact into content. The Marseille
calibrated artifact remains one section, with the raw-closure gap stated and
partially decomposed (frame-convention fix, chi-sign fix, remaining physics
gap under active investigation).

## Skeleton (framing B; sections collapse gracefully into framing A if chosen)

1. **Introduction** — polarized twilight RT: why backward MC; why validation
   of polarization chains is hard (rotation-convention bugs are invisible to
   unpolarized single-scatter tests); prior twilight models.
2. **Model** — geometry (spherical shells, SZA>90 tangent paths), Stokes
   transport, deterministic single + second-order control variates, MC
   higher orders with twilight-guided branching and the recursive branch cap,
   estimator options (mean vs median-of-means, robust groups), knobs with
   frozen-reproducibility defaults. Source: `monte_carlo_cpp/README.md`,
   `MODEL_ASSUMPTIONS.md`.
3. **Validation suite and reconciliation protocol** — the seven claim
   families; per-case fresh-load gate mechanics; RNG scheme (case-id
   independent, bit-reproducible reruns); the reconciliation standard
   (`tools/reconcile_validation_numbers.py`, 0/7 families missing
   2026-08-03). Table 1 = the final scoreboard: convergence PASS under
   declared cap=1 deviation / DISORT exact / Zawada single+multiple ~1.000 /
   IPRT legacy-fail vindicated under fix / Rozenberg FAIL (under
   investigation) / Koomen mixed.
4. **The chi-sign regression: archaeology and fix** — the centerpiece.
   eventMuellerMatrix applied R(-chi_in) with R(+chi_out) though both chi are
   measured about backward ray axes; invisible to unpolarized inputs (single
   scatter exact), corrupts any polarized input at an event. Independent
   dipole double-scatter ground truth (45 geometries, exact to 1e-9);
   fix-on IPRT metrics equal the historical claims to 4+ sig figs
   (median_dolp_abs = 0.0055043 exactly); committed history always carried
   the regression (bit-identical from-source rebuild at the claim commit);
   Zawada all-orders improves ~3x as an untouched external cross-check.
   Artifacts: `*_chifix*` in `paper/tables/reproducibility/`.
5. **Non-terminating recorded configurations and declared deviations** —
   uncapped twilight branching diverges (branch factor 2 per order-3+ event);
   the recorded default gate ran 9 days without completing one stage while
   the cap=1 gate finished in 21 s on the same binary; historical
   convergence/Rozenberg/Koomen quotes are unproducible from recorded
   configs; the declared-deviation protocol (cap=1 in-file with provenance
   comments) as the honest reconciliation device.
6. **The open twilight gap** — raw Marseille signed DoLP +0.14 legacy /
   +0.19 corrected; Rozenberg normalized profile falls off the horizon peak
   3-15x faster than the reference with single-scatter dominance where the
   real sky is multiple-scatter dominated; guard ladder shows order content
   saturated by ~order 8-16 (over-polarization formed at low orders, not a
   truncation artifact); King-factor molecular depolarization quantified
   (knob `rayleigh_depolarization_factor`, DoLP(90) = (1-rho)/(1+rho));
   candidate mechanisms under test (guided-estimator MS yield via no-guiding
   control; missing refraction; aerosol profile). RESULTS PENDING: jobs
   787251/787252 + build_fix2 leg.
7. **Calibrated Marseille full-field artifact** — the frozen April result,
   scoped exactly as CLAIMS.md dictates (calibrated pipeline closure, not raw
   closure); the raw-vs-calibrated distinction becomes a virtue in framing B.
8. **Reproducibility statement** — frozen tiers, knob-gated changes with
   bit-identical defaults, archived reports + reconciliation outputs, Slurm
   job provenance.

## Figure inventory

Existing (frozen): Marseille calibrated main panel (2026-04-30).
To generate (sources on disk):
- F1 chi-fix mechanism diagram + dipole ground-truth error table (45
  geometries; scratch probe archived in the progress log).
- F2 IPRT A1 DoLP: reference vs legacy vs fix (fix-on = claims to 4+ figs);
  CSVs in `paper/tables/reproducibility/`.
- F3 order-resolved twilight decomposition: guard ladder + first/second/
  higher fractions (comparison CSVs, `*_guard{8,16,32}`).
- F4 Rozenberg normalized zenith profiles, model vs reference, annotated
  with first/second/higher fractions (`measurement_rozenberg_hminus6_comparison.csv`).
- F5 Marseille signed-DoLP sky map, legacy vs corrected chain (full-field
  chifix run, job 787251, landing).
- F6 the 9-day/21-second control pair (stdout timelines; sacct records).

## Claims discipline

`paper/CLAIMS.md` remains binding for the calibrated artifact. New claims to
draft (framing B): chi-fix vindication claim (exact numbers, no adjectives),
declared-deviation convergence claim, open-gap claim with signed values.
Forbidden: presenting capped/deviation numbers as the recorded-config values;
presenting the calibrated Marseille pass as raw closure.

## Blocking inputs before full drafting

1. Framing A vs B.
2. Depolarization verdict (jobs in flight: low-guard + King-factor leg,
   no-guiding Rozenberg control 787252, full-field corrected chain 787251).
3. Chi-fix adoption decision + the new convergence-under-fix anomaly
   (candidate gate 787253: corrected chain fails the default gate's
   convergence family at 256 photons - heavier-tailed higher-order intensity;
   needs diagnosis before any production freeze).
