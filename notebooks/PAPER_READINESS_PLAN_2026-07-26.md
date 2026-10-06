# Paper Readiness Plan

Prepared: July 26, 2026
**Status (updated 2026-07-26 evening): implementation underway — see [PLAN_IMPLEMENTATION_PROGRESS_2026-07-26.md](PLAN_IMPLEMENTATION_PROGRESS_2026-07-26.md).**
A1 done (5/5 ctest pass); A2/A3 running in background; C1 done (median-of-means suppresses higher order ~95%, all 683 directions); B1 done (intensity shape error is ⅓ estimator bias; gate never compares absolute scale); B4 done (raw AoP misfit is a model-side polarization-frame reflection: 43°→14.4° median after mirror+const — fix-site named). B5 partially satisfied by exact offline re-scoring of the full frozen field.

Based on: [REALISTIC_LIGHT_SCATTERING_READINESS_AUDIT.md](audits/REALISTIC_LIGHT_SCATTERING_READINESS_AUDIT.md) (overall grade 70/100, C-; verdict: not yet ready, substantial technical/validation work remains).

## Bottom line

The audit found two things worth separating clearly, because they call for different kinds of work:

1. **A confirmed-broken evidentiary claim**: the flagship Marseille "calibrated" result is a mathematical tautology (per-point calibration with enough free parameters to force exact agreement, regardless of physics). This isn't a validation gap to close — it's a claim to stop making, and replace with something real.
2. **A confirmed-large, undiagnosed physics gap**: the raw (uncalibrated) model misses the same Marseille reference by ~9 orders of magnitude in intensity and ~44° in AoP. This is genuinely open science, not a checklist item.

Everything below is sequenced around one rule learned the hard way by this project already: **diagnose cheap, confirm expensive.** The strict full 683-direction Marseille run has, historically, cost 72+ hours per higher-order sample and produced zero durable progress in some attempts (`notebooks/MARSEILLE_JOB_FAILURE_REPORT_2026-04-17.md`). Every diagnostic phase below uses the existing cheap subset configs (`frozen_marseille_twilight_20220815_191413z_measurement_strict_subset.cfg`, `..._profile_subset.cfg`, `..._tiny_subset.cfg`) or reuses already-computed artifacts. Full-field strict reruns happen exactly twice in this whole plan: once to confirm a fix at the end of Phase B, and once for the final paper-grade number.

---

## Guiding principles

- **Re-establish ground truth before diagnosing anything.** We don't currently know if the general benchmark suite (Zawada, DISORT, IPRT) still passes, because its result artifacts don't exist on disk. Find out first — it's cheap and it tells you how much of the "the core machinery works" story you get to keep.
- **Separate "pipeline plumbing works" from "physics is correct."** The calibrated-Marseille work already done is legitimate evidence of the former. Don't let it get recharacterized as the latter in any future doc or manuscript text.
- **The paper direction is a decision gate, not a starting assumption.** Section E below lays out three candidate papers with explicit trigger conditions. Don't pick one before the evidence is in.
- **Every phase ends with a stated success criterion**, not just a task list, so "done" is checkable.

---

## Phase A — Re-establish the evidentiary baseline

**Targets:** audit Finding M-2 (general validation suite artifacts missing from disk and git).
**Why first:** cheap, fast, and everything downstream implicitly assumes the answer.

| # | Task | Notes |
|---|---|---|
| A1 | Rebuild with `build_cluster_gcc8` and rerun `ctest` | Confirms the 5 unit-style physics tests are still green (they were as of 2026-05-01; should take seconds) |
| A2 | Run `ValidationRunner` against `default_clear_sky.cfg` | Regenerates `results/validation/validation_report.json`. Not the expensive strict path — should be tractable interactively or as a short batch job |
| A3 | Rerun both Zawada spherical-vector benchmark configs, DISORT scalar, IPRT A1 vector, Rozenberg, Koomen | These are the strongest non-Marseille evidence the project claims to have; right now that claim is unconfirmable |
| A4 | **Durably preserve the outputs.** Don't let them fall back into the gitignored `results/` void. | Recommend copying final JSON/txt summaries into `paper/tables/reproducibility/`, mirroring the pattern already used for the Marseille artifacts, or carve a narrow `.gitignore` exception for validation summary files specifically |
| A5 | Reconcile against the numbers currently quoted in `STATUS_README.md` / `MODEL_ASSUMPTIONS.md` / `monte_carlo_cpp/README.md` | Update any that drifted; note any that didn't reproduce |

**Success criterion:** `validation_report.json` and both Zawada smoke outputs exist, are committed or otherwise durably preserved, and their numbers are either reconfirmed or honestly updated.
**Effort:** small (hours–1 day). **Blocking:** not strictly, but do it first anyway — it's the cheapest phase and de-risks everything after it.

---

## Phase B — Diagnose the raw Marseille mismatch

**Targets:** audit Finding M-1 — the real bottleneck of this entire plan.
**Ground rule:** everything here runs on the strict-subset or profile-subset configs. No full 683-direction reruns until B5.

| # | Task | Notes |
|---|---|---|
| B1 | Units/normalization pass | The reference intensity looks like near-raw camera DN (~10²–10³); the model outputs physical spectral radiance (~10⁻⁷). Check whether `instrument_response.csv` / the IMX250MYR proxy is actually applied anywhere as a DN conversion, or only as a spectral weight in `buildBands`. Back out the single best-fit *constant* conversion factor, then look at what's left — that residual is the real, un-explained-by-units mismatch, since the current per-point gain (which varies ~50× across the sky) conflates the two |
| B2 | Refraction test | The model has zero atmospheric refraction, in exactly the below-horizon geometry where refraction is best known to matter. Prototype a simple refraction correction to the solar-geometry/line-of-sight calculation and rerun the strict subset to see whether raw DoLP/AoP error drops |
| B3 | Aerosol/atmosphere input sensitivity | The frozen case's aerosol optics use `aeronet_inversion_usage_mode = same_day_fallback`, not a strict within-window match. Bound how much this could plausibly explain by perturbing/swapping the aerosol input and re-scoring the subset |
| B4 | Marseille reducer geometry audit | Independently re-derive the `rotation.npy`-based Q/U basis rotation in `reduce_marseille_case.py` from `dataset_readme.md`'s own description, rather than trusting the current implementation. Test against a synthetic case with a known polarization angle |
| B5 | Structured ablation + one confirmation run | Once B1–B4 produce candidates, toggle each independently (then combined) on the strict subset and track normalized RMSE / DoLP / AoP error. Once a fix (or fix-combination) is identified, run it once, full-field, strict, to confirm |

**Success criterion:** raw (uncalibrated) comparison shows a materially reduced, *understood* mismatch — ideally down near the scale of the already-documented Gal-case failure (a few percent to ~0.05–0.1 DoLP-scale error) or better — with a stated explanation for what closed the gap.
**Effort:** medium–large; open-ended diagnostic science, realistically 1–4+ weeks depending on what's found. This is the actual critical path of the whole plan.
**Decision point (end of Phase B):** gap closes substantially → proceed toward the physics-validation paper (E-1). Gap doesn't close → fall back to the methods/software paper (E-2), with Marseille reported as an honest open limitation, not a validation case.

---

## Phase C — Numerical estimator validation

**Targets:** Finding M-3 (unvalidated median-of-means higher-order estimator) and the missing strict-settings convergence check.
**Can run in parallel with Phase B** — it doesn't depend on B's outcome.

| # | Task | Notes |
|---|---|---|
| C1 | Median-of-means vs. plain mean, on a handful of representative directions | Reuse the profile-subset directions (near-zenith, mid-zenith, horizon-skimming) — their per-direction costs are already characterized in `profile_subset_runtime_probe.log` |
| C2 | If they diverge materially, fix it | Either raise samples-per-group or swap the estimator for production use |
| C3 | Run the existing convergence gate (half- vs. full-photon-budget) at the *actual strict Marseille configuration*, not just `default_clear_sky.cfg` | The checked-in gate currently only exercises the default config |

**Success criterion:** documented evidence the higher-order point estimate isn't materially biased, plus a convergence check that actually exercises strict Marseille settings.
**Effort:** small–medium (reuses existing subset infrastructure).

---

## Phase D — Holdout validation of the calibration (conditional)

**Targets:** Finding C-1. **Only pursue this if Phase B narrows but doesn't fully close the raw gap** — if B fails outright, skip to E-2; if B fully closes the raw gap, this phase becomes unnecessary (a genuinely-closing raw model is strictly stronger and simpler than any calibrated one).

| # | Task | Notes |
|---|---|---|
| D1 | Split the 683 directions into fit/holdout sets | Use spatial blocks (e.g. the existing `region_*` groupings already computed in the region-summary output), not random points — neighboring sky bins are correlated |
| D2 | Fit whatever calibration (or replacement) survives Phase B on the fit set only | |
| D3 | Score the held-out set and report that number, honestly, as the real result | This is the number that would actually mean something — versus the current one, which is true by construction |

**Success criterion:** a calibration (if kept at all) with demonstrated out-of-sample skill, reported instead of the current in-sample-by-construction number.
**Effort:** small — this is pure data analysis on artifacts that already exist (`measurement_model_quality_calibration.csv`, the comparison CSV), no reruns needed.

---

## Phase E — Decide the paper direction (explicit gate)

Don't resolve this before Phase B reports. Three candidates, each with a trigger condition:

| Direction | Trigger | Notes |
|---|---|---|
| **E-1: Twilight polarization validation paper** | Phase B closes the raw gap to a defensible level | Strongest, most interesting outcome — an actual reproduction/validation result against real measurements |
| **E-2: Methods/software paper** | Phase A confirms the general benchmarks (Zawada, DISORT, IPRT) genuinely pass — independent of how Phase B resolves | The safe minimum-viable deliverable; can be written even if Marseille never fully closes, with that framed explicitly as a limitation / future work |
| **E-3: Parameter-sensitivity or physical-interpretation study** | Only after E-1 or E-2's groundwork | Needs net-new sweep campaigns beyond anything that currently exists (only one null-result surface-sensitivity check exists today). Treat as a second-paper goal, not the first target |

**Recommendation:** structure Phases A–D so that E-2 is fully supported regardless of how Phase B turns out — it's the floor. E-1 is the upside case.

---

## Phase F — Manuscript production

Only after E is decided.

| # | Task |
|---|---|
| F1 | Uncertainty quantification for whatever quantitative claims the chosen paper makes (feeds from Phase C) |
| F2 | Novelty positioning against the still-unpublished Poughon et al. dataset paper — periodically check whether the "BMC Research Notes" publication has appeared, since that changes what's citable and how this work should be framed relative to it |
| F3 | Regenerate `paper/figures/marseille_calibrated/*` once the calibration story changes post-Phase B/D, via the existing `paper/reproduce.py --regenerate-plots` |
| F4 | Update `paper/CLAIMS.md` / `text_snippets.md` — only remove the "not independent raw first-principles closure" hedge if Phase B/D evidence actually earns it |
| F5 | Standard figure polish, drafting, internal review, submission |

---

## At-a-glance roadmap

| Phase | Goal | Depends on | Rough effort | Blocking for paper? |
|---|---|---|---|---|
| A | Confirm general benchmark suite still passes | — | Small (~1 day) | No, but do it first |
| B | Diagnose raw Marseille mismatch | A (informative, not blocking) | Medium–large (1–4+ wks) | **Yes — the critical path** |
| C | Validate higher-order estimator | — (parallel with B) | Small–medium | Yes, for any quantitative claim |
| D | Holdout-validate calibration | B (conditional) | Small | Only if B partially succeeds |
| E | Choose paper direction | A, B, C, D | — (decision, not work) | Yes |
| F | Write the paper | E | Varies by direction | — |

---

## What not to do

- **Don't rerun the full 683-direction strict Marseille case as a diagnostic tool.** That's the exact trap this project already fell into (see `MARSEILLE_JOB_FAILURE_REPORT_2026-04-17.md`). Full-field strict runs happen twice, total, in this plan: once to confirm a hypothesized fix (end of Phase B), once for the final paper number.
- **Don't report calibrated Marseille numbers as physics validation in any manuscript text** until Phase D's holdout numbers exist to back that framing up.
- **Don't re-enable the ad hoc twilight tuning knobs** (`boostUnpolarizedIntensity`, `twilight_order_depolarization`, etc.) as a shortcut during Phase B. If Phase B finds a real fix, these should be deleted from the codebase rather than left dormant as a temptation for the next deadline crunch.
- **Don't let housekeeping block the critical path.** The stale Windows-path links in `STATUS_README.md`/`MODEL_ASSUMPTIONS.md` and the coarse git commit messages (audit Findings M-min-1, M-min-2) are real but low-stakes — fix them opportunistically alongside any other phase, not as a gating task.
