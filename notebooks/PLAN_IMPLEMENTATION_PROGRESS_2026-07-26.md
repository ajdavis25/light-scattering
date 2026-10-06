# Paper Readiness Plan — Implementation Progress

Prepared: July 26, 2026 (evening session)
Plan: [PAPER_READINESS_PLAN_2026-07-26.md](PAPER_READINESS_PLAN_2026-07-26.md)
Audit basis: [audits/REALISTIC_LIGHT_SCATTERING_READINESS_AUDIT.md](audits/REALISTIC_LIGHT_SCATTERING_READINESS_AUDIT.md)

## Headline results of this session

The plan's diagnostic phases produced three major findings, each closing or reframing one of the audit's open problems. **All three raw-Marseille misfit axes now have distinct, quantified, named causes:**

| Misfit axis (raw, audit M-1) | Cause found | Evidence | Fix class |
|---|---|---|---|
| AoP: median 43° error, near-random vs reference | **Polarization-frame reflection (handedness/sign convention) on the model side of the comparison** — not physics | Mirror + one constant rotation (+82.8°) collapses median AoP error 43.0° → **14.4°** (6.8° near zenith, 11.9° in DoLP>0.15 bins); mirrored model locks to independent single-scatter geometry at **0.883** | Convention fix in the measurement-comparison path (post-processing-level; no new solver runs required to validate) |
| Intensity shape: 0.188 masked normalized RMSE | **~1/3 is median-of-means estimator bias** (audit M-3, now proven); remainder has strong scattering-angle structure | Re-scoring the frozen field with the unbiased plain-mean higher-order component: masked RMSE 0.188 → **0.122**; residual correlates with scattering angle (r = −0.53: model relatively too bright sunward) | Estimator change + then genuine physics/inputs work on the sunward-brightness residual |
| DoLP: +0.148 mean over-polarization | **Genuine physics deficit** — model is essentially an over-polarized Rayleigh sky; NOT estimator, NOT convention | Estimator swap barely moves it (median 0.187 → 0.183); model DoLP correlates 0.87 with pure Rayleigh single-scatter shape vs 0.58 for the measurement; model DoLP reaches 0.93 vs measured max 0.58 | Real modeling work: depolarization pathways (aerosol depol, multiple scattering content, possibly refraction geometry) |

A fourth structural clarification: **the audit's "9 orders of magnitude" intensity headline needs reframing.** The measurement gate peak-normalizes each field independently (`measurement_case_main.cpp:1635-1638`), so absolute scale was never compared — the comparison is DN-vs-radiance with no conversion anywhere, by design invisible to the gate. The scientifically meaningful intensity errors are the shape numbers above. The log-log slope of reference vs model is 0.51 (the model has ~2× the measurement's dynamic range in log space), consistent with suppressed multiple scattering; a slope-0.5 law explains ~10× more variance than any constant gain (R² 0.43 vs 0.04).

---

## Phase-by-phase status

### Phase A — Re-establish the evidentiary baseline: IN PROGRESS (launched, partially complete)

- **A1 DONE.** `ctest` rerun fresh: **5/5 pass** (phasefunctions, surface, polarization, atmosphere, wavelength; 0.29 s total; log in `monte_carlo_cpp/build_current/Testing/Temporary/LastTest.log`, dated this session). The stale `LastTestsFailed.log` confusion from the audit is superseded.
- **A2/A3 LAUNCHED, RUNNING.** New sequential driver `monte_carlo_cpp/tools/run_phase_a_validation.sh` (16 threads, `nice`, from repo root) runs: (1) the default gate — convergence + DISORT scalar + IPRT A1 vector + Rozenberg + Koomen; (2) a **new** standalone smoke wrapper for the Zawada single-scatter benchmark (`config/benchmark_zawada_spherical_vector_single_smoke.cfg`, created this session — none existed); (3) the existing Zawada all-orders smoke wrapper. The driver archives `validation_report.json` and runner stdout under per-case names after each run, fixing the clobber hazard. Driver log: `monte_carlo_cpp/results/validation_phase_a_driver.log`; outputs land in `monte_carlo_cpp/results/validation/`.
- **A4/A5 PENDING** the runs completing: copy summaries into `paper/tables/reproducibility/` and reconcile against the numbers quoted in `STATUS_README.md` / `MODEL_ASSUMPTIONS.md` / `monte_carlo_cpp/README.md`.
- **New gap found during A:** the benchmark gate loads signed Q/U from the vector references (Zawada, IPRT/MYSTIC both carry signed `q,u` columns) but **compares only normalized intensity and |DoLP|** (`ValidationQA.cpp`, benchmark path). Signed Q/U — the exact quantity implicated in the AoP convention finding below — has never been scored by any benchmark. Adding signed-Q/U benchmark metrics is now the decisive convention-anchoring test and a permanent regression guard. (Follow-up task; small code change in `ValidationQA.cpp`.)

### Phase C1 — Estimator comparison: DONE (and decisive)

New durable tool: `monte_carlo_cpp/tools/analyze_higher_order_estimator.py`. It exploits the fact that the frozen robust_r2 run's 683 per-direction checkpoint files store **both** estimators' ingredients (plain running mean `higher_mean_*` + all 16 robust-group sums). Reconstruction was verified exact against the frozen comparison CSV (0/683 mismatches, worst relative deviation 4.7e-11), so these are the frozen field's actual numbers, not approximations:

- Median-of-means sits **below** the plain mean in **683/683 directions**; median suppression of the higher-order intensity component ≈ **95%**.
- Field-summed intensity: median-of-means field is **21.4% dimmer** than the plain-mean field (ratio 0.786). Per-direction total-intensity difference: median 6.8%, p95 53%, max 97%.
- Estimator choice alone shifts DoLP by up to 0.19 (p95) and AoP by up to 28° (p95) in tail directions.
- The plain mean is not trustworthy either at this budget: from the stored Welford `m2` accumulators, its relative standard error on the higher-order component is **median 62%** (501/683 directions above 50%) at 256 samples.
- **Conclusion (audit M-3 resolved):** at the frozen run's budget, *neither* estimator reliably measures the higher-order component; the frozen field carries a systematic negative multiple-scattering bias. Because suppressing multiply-scattered (depolarizing) light over-polarizes the sky, this bias pushes in the same direction as the dominant raw DoLP failure — though the estimator-swap re-score shows it accounts for only a small part of the DoLP misfit (it accounts for a large part of the intensity-shape misfit instead).
- Per-direction output: `monte_carlo_cpp/results/measurement_case_reports/higher_order_estimator_comparison.csv`.

### Phase B1 — Units/normalization decomposition: DONE

New durable tool: `monte_carlo_cpp/tools/analyze_raw_marseille_misfit.py`; per-row output `monte_carlo_cpp/results/measurement_case_reports/raw_misfit_decomposition.csv`. Key numbers beyond the headline table above:

- Best single constant gain 2.52e9 (10^9.40); residual after one constant: 0.249 dex (×1.78 at 1σ) — so no single unit constant explains the spread, but see the slope result: free log-log fit slope **0.514**, R² 0.43 (vs 0.04 const-only).
- Reference "intensity" values run 339–23,181 with p95 = 4,136 — consistent with the dataset's documented 12-bit DN bit-shifted encoding, not physical radiance; no DN→radiance conversion exists anywhere in the comparison chain (and none is needed for the peak-normalized gate).
- Masked (≥5% of peak) shape residual: model relatively too bright at small scattering angles (r = −0.525 vs scattering angle); 287/323 masked bins have normalized model > normalized reference.
- DoLP misfit worst near zenith (median |ΔDoLP| 0.36 for zenith < 30°) — the region where measured twilight DoLP should be most multiple-scattering-diluted.

### Phase B4 — Reducer/convention audit: DONE (diagnosis complete, fix-site named)

- The reducer's own math is clean: `rotate_linear_stokes` (reduce_marseille_case.py:81-91) is a proper det=+1 basis rotation (verified analytically: `aop' = aop − angle`, no reflection), and the reduced reference field locks to independent geometric single-scatter expectation at 0.589 with the *same* handedness.
- The model, as scored, locks to the reference at only **0.030** (uncorrelated); the **mirrored** model locks at **0.573** to the reference and **0.883** to geometric expectation. The full relation is `reference ≈ 82.78° − model_aop` — a single reflection, equivalently a sign flip of model Q in the comparison basis (plus a small −7° residual constant, within the post-fix residual).
- **Fix site:** the measurement scorer writes `model_q/model_u` directly from the solver's outward-ray (backward-tracing) basis with no comparison-basis conversion (`measurement_case_main.cpp`; duplicated in `ValidationQA.cpp` measurement path), while the reduced measurement is quoted in an incoming-light convention. One reflection at that interface reconciles them.
- This retroactively explains the long-chased "AoP flip sectors": under a reflection about the ~41° axis, the documented 45°/215° relative-azimuth sectors (median ≈ 16.7° in the interactive run) are simply the sky regions nearest the reflection axis. It also explains why the earlier "azimuth-convention sweep" found nothing: it rotated the field's *position*, which cannot detect an angle-*value* reflection.
- Post-fix residual is zenith-structured (median 6.8° near zenith → 43° for zenith 60–90°), consistent with the C1 multiple-scattering deficit dominating near the horizon, where AoP is also least well-defined (low DoLP).

### Phase B5 (candidate-fix confirmation): PARTIALLY SATISFIED OFFLINE

The convention fix and the estimator swap were both validated by **exact offline re-scoring of the full frozen 683-direction field** (stronger than the planned strict-subset A/B, since it uses the complete real run at zero compute cost). What offline re-scoring cannot do: produce a new, higher-budget, convention-fixed solver run. That remains the one genuinely open compute task (see next steps).

### Phases C2/C3, D, E, F: NOT STARTED (correctly blocked)

- C2 (estimator fix in production) and the signed-Q/U benchmark metric are now well-specified small code changes; C3 (convergence at strict settings) should follow the estimator decision, since convergence of a biased estimator is moot.
- D (holdout calibration) — the diagnosis suggests the row-wise calibration may become unnecessary: with the convention fix, AoP drops to ~14° median raw; the remaining DoLP/shape work is physics, not calibration.
- E (paper direction): evidence now favors an upgraded outcome — see below.

---

## Updated read on the paper direction (Phase E, preliminary)

Today's findings move the project meaningfully toward the **E-1 (validation-paper)** trigger rather than the E-2 fallback:

- The "near-random AoP" — the scariest audit finding — is a comparison-convention artifact, not model failure. The model's polarization geometry is *good* (0.883 lock with theory after handedness normalization).
- The intensity shape error is one-third estimator artifact, with a physically interpretable sunward residual.
- The genuinely open physics is now narrow and nameable: **insufficient depolarization/multiple-scattering content** (DoLP +0.12–0.15 over-polarization surviving all corrections) — with the C1 sampling deficit (62% median SE on higher order) as the immediate suspect to eliminate first, since the higher-order component is currently both biased *and* noisy.

None of this is paper-ready until: (1) Phase A confirms the benchmark suite (running), (2) the convention fix + estimator fix land in code with signed-Q/U benchmark guards, and (3) one convention-fixed, adequately-sampled strict run (or batch) replaces the frozen artifact.

## Immediate next steps (in order)

1. **(Running)** Let the Phase A driver finish; then A4/A5: preserve artifacts to `paper/tables/reproducibility/`, reconcile doc numbers.
2. **Code: DONE (this session, later).** All three changes implemented, config-gated with defaults preserving current behavior so nothing changes silently:
   - `higher_order_estimator = median_of_means | mean` (`MonteCarloDriver.hpp/.cpp`; selection at the `medianOfMeans` call site; value-validated at config load; echoed in the result JSON).
   - `measurement_model_polarization_convention = legacy_outward_ray | incoming_light` (`MonteCarloDriver.hpp/.cpp` for config plumbing; `applyIncomingLightConvention` = Q→−Q applied once, immediately post-solve, in **both** scorer paths — the direct full-solve path and the batched single-direction path in `measurement_case_main.cpp`, plus `ValidationQA.cpp`'s measurement path; echoed in the measurement report txt for provenance).
   - Signed-Q/U diagnostic metrics in the benchmark path (`ValidationQA.cpp`): `diag_signed_qn_median_abs`, `diag_signed_un_median_abs`, `diag_q_sign_agreement_frac`, `diag_u_sign_agreement_frac` — never gating (pass always true), computed from signed q/I, u/I against the vector references. These anchor the solver's Stokes convention against MYSTIC/Zawada and guard the convention permanently.
   - New configs (frozen originals untouched): `..._measurement_tiny_smoke_conventionfix.cfg` (A/B verification) and `..._measurement_strict_subset_fixed.cfg` (go-forward strict-subset confirmation, convention fix + mean estimator).
   - Built in isolated `monte_carlo_cpp/build_fix/` (leaves `build_current` untouched for the running Phase A jobs); **ctest 5/5 pass on the modified build**.
3. **Verify: transform invariants CONFIRMED end-to-end (2026-07-27).** The original login-node verification attempt was killed by a session restart after ~5 h (lesson recorded below), and was replaced by a purpose-built micro A/B pair — 3 directions, `single_scatter_only=true`, since the convention transform is a deterministic invariant that does not need scattering-order realism. Result (`tools/verify_convention_ab_pair.py` logic run on `..._measurement_micro{,_conventionfix}` comparison CSVs, both produced by the real `build_fix` binary in <0.2 s each): **all invariants PASS exactly** — `model_q` flips sign to the last bit, `model_u`/intensity/DoLP unchanged, `model_aop = 90° − legacy_aop`, and both report files echo their `measurement_model_polarization_convention` value correctly.
4. **Compute moved to Slurm (2026-07-27), replacing the failed login-node runs.** Two jobs submitted (**outcome 2026-07-29: both hit 48 h walltime — see the dated section at the end of this file for the post-mortem, the replacement jobs 779980–779983, and the IPRT finding**):
   - **778296** `phase_a_validation` (compute1, 16 cpus, 48 h): regenerates the validation artifacts with the **unmodified** solver (`build_cluster_gcc8`), fastest-first (Zawada single → Zawada multiple → default gate), per-stage clobber-safe archiving. Reconciliation against historically claimed numbers is one command when it lands: `python3 monte_carlo_cpp/tools/reconcile_validation_numbers.py`.
   - **778297** `marseille_strict_subset_fixed` (compute1, 16 workers, 48 h): the Phase B5 confirmation run — 48 strict-subset directions under the **fixed** build (`build_fix`: incoming-light convention + mean estimator), via the resumable per-direction checkpointing batch runner (`--higher-order-block-size 1 --resume`), so walltime kills lose at most one higher-order block. Offline predictions to beat: raw median AoP ≈ 14° (was 43°), masked shape RMSE ≈ 0.12 (was 0.19).
   - Slurm scripts: `monte_carlo_cpp/slurm/phase_a_validation.slurm`, `monte_carlo_cpp/slurm/strict_subset_fixed.slurm`.
5. **Then** revisit Phase D/E with the new raw numbers.

### Operational lessons from the failed first attempt (2026-07-26 evening)

- Login-node background runs do not survive Claude session restarts; anything longer than ~minutes belongs in sbatch (this is also just correct multi-user etiquette). Both original runs (default gate at 4.5 h, tiny A/B pair at ~6 h) died with zero recoverable output.
- **New runtime intel from the killed tiny run's log:** at only 2 photons/bin and a 430–490 nm band, one tiny-subset direction spent **4.2 hours** in the higher-order stage while most others took seconds. The twilight higher-order branching cost remains extremely direction-dependent — any future budget planning (including the C1-motivated sample-count increase) must budget per-direction, not per-case, and the per-direction checkpointing path is the only safe execution surface for it.

## Artifacts created this session

| Path | What |
|---|---|
| `monte_carlo_cpp/tools/run_phase_a_validation.sh` | Sequential Phase A validation driver (clobber-safe archiving) |
| `monte_carlo_cpp/config/benchmark_zawada_spherical_vector_single_smoke.cfg` | New standalone smoke wrapper for the Zawada single-scatter benchmark |
| `monte_carlo_cpp/tools/analyze_higher_order_estimator.py` | C1 estimator comparison (exact, checkpoint-based) |
| `monte_carlo_cpp/tools/analyze_raw_marseille_misfit.py` | B1 raw-misfit decomposition |
| `monte_carlo_cpp/results/measurement_case_reports/higher_order_estimator_comparison.csv` | Per-direction estimator deltas |
| `monte_carlo_cpp/results/measurement_case_reports/raw_misfit_decomposition.csv` | Per-direction misfit decomposition (gain, shape, DoLP, AoP, scattering angle) |
| `monte_carlo_cpp/results/validation/` (populating) | Regenerated validation artifacts (Phase A) |
| This file | Session record |

No frozen artifacts were modified; all new outputs are additive.

---

# Session update 2026-07-29: job post-mortems, replacement jobs, and the IPRT finding

## Post-mortem of jobs 778296 / 778297 (both TIMEOUT at 48 h, by design partially productive)

**778296 `phase_a_validation`:** the two Zawada smokes landed in the first 35 s and reconcile cleanly
(see below). The default gate then ran silently for 47.9 h: `runValidationSuite` executes the
convergence stage FIRST — two full-sky runs of `default_clear_sky.cfg` (coarse = photons/2, then
fine), 684 bins x 46 wavelengths, all orders at solar zenith ~97 deg — before any benchmark case is
evaluated. The stdout silence is a buffering artifact (all 19 "Evaluating zenith row" lines print
up front while building the direction list; the sky work happens afterwards inside the
OpenMP-parallel `sampleSkyDirections`, which prints nothing without a callback). Consequence:
none of DISORT/IPRT/Rozenberg/Koomen were ever reached.

**778297 `marseille_strict_subset_fixed`: the strict-subset tier is computationally infeasible, full stop.**
0/48 directions completed; all 16 workers wrote their initial checkpoints in the first minute and
then sat inside their first higher-order sample for 48 h (only 4/16 even logged band 1/25).
Root cause found by diffing the config tiers: the frozen production run
(`__full_branchcap_robust_r2`) sets `higher_order_recursive_branch_cap=1`; the `strict_subset`
tier REMOVES the cap, so 3-way twilight branching (`twilight_higher_order_branches=3`) recurses
uncapped under `max_events_guard=64` — combinatorial path explosion at 96 deg solar zenith.
This also retroactively explains the killed tiny-subset run's 4.2 h/direction observation.
**The branchcap tier is the production tier; it is also the tier every offline prediction was
computed from, so it is the correct confirmation vehicle anyway.**

## Replacement jobs (submitted 2026-07-29, monitor armed)

| Job | Script | Partition | What |
|---|---|---|---|
| 779980 | `slurm/phase_a_split_cases.slurm` | compute1, 36 cpu, 3 d | split wrappers: DISORT -> IPRT -> Rozenberg, sequential, per-stage archiving |
| 779981 | `slurm/phase_a_split_koomen.slurm` | compute2, 36 cpu, 5 d | Koomen split wrapper (heaviest case, own walltime) |
| 779982 | `slurm/validation_default_gate.slurm` | compute2, 36 cpu, 9 d | full unmodified default gate (the real convergence artifact; ~6 d by core-hour extrapolation) |
| 779983 | `slurm/branchcap_subset_fixed.slurm` | compute2, 32 workers, 6 d | Phase B5 confirmation at the FROZEN branchcap budget + the two fix knobs, `--resume`-safe |

Split-wrapper design (`config/validation_split_{disort,iprt,rozenberg,koomen}.cfg`): clone of
`default_clear_sky.cfg` with a 1x1-bin, single-wavelength top-level sky (the convergence prefix
becomes trivial and meaningless — convergence claims may ONLY be reconciled from the full gate)
and exactly one `benchmark_case_config`/`measurement_case_config` entry. This is rigorous because
`evaluateBenchmarkConfig` loads the case config FRESH (`loadSimulationConfig(configPath)`) — no
inheritance from the wrapper. Verified empirically: wrapper runs reproduce case metrics bit-for-bit
across builds and thread counts.

**Subset-specific predictions for 779983** (recomputed within the 48-direction field —
full-field numbers do NOT transfer; recorded also in the slurm script header):
- raw median |AoP err| 52.5 deg -> **21.2 deg** with `incoming_light` (subset oversamples hard regions vs full-field 43 -> 14.4)
- masked (31/48 rows) peak-normalized shape RMSE: MoM 0.2531 -> mean 0.2554 — **the estimator swap does NOT improve shape on this subset** (the 0.188 -> 0.122 gain was a full-field property). Do not use shape RMSE as the estimator pass criterion here; use per-direction consistency with frozen totals re-aggregated as mean.
- DoLP signed err median **+0.134** persists (genuine physics gap).

## Phase A reconciliation status (fresh vs claimed; `tools/reconcile_validation_numbers.py`, now label-aware)

| Family | Status | Fresh vs claimed |
|---|---|---|
| DISORT scalar | **REPRODUCED (5 significant figures)** | 0.0188114/0.0250572 vs 0.0188119/0.0250578 |
| Zawada single | reproduced (MC-noise level) | ratios 1.00–1.19 |
| Zawada multiple | reproduced | ratios 0.99–1.00 |
| IPRT A1 vector | **NOT REPRODUCED — see finding below** | intensity 3.9x claim; DoLP 16.8x claim, FAILS gate |
| Rozenberg | **NOT REPRODUCED** (780007, capped rerun; see 19:45Z note) | RMSE 0.264325 vs 0.0849321 claimed — 3.1x, FAILS gate |
| Koomen | **MIXED** (779981, capped, 19 s; see 20:45Z note) | RMSE 0.102507 vs 0.0470789 claimed (2.2x, but still passes gate); p95 DoLP 0.0309 vs 0.0367 claimed (fresh BETTER) |
| convergence | pending (779982, only valid source) | — |

## FINDING (major): IPRT A1 vector DoLP — claim unsubstantiated at its own commit; real multiple-scatter polarization defect

Fresh numbers (identical to the last digit across `build_current`, `build_cluster_gcc8`,
`build_fix`, AND a from-source build of commit `232f67f` — the very commit that introduced the
claim text, the case config, and the reference data — at 8/9/16/36 threads):
`median_intensity_rel = 0.0220254` (claimed 0.00561889), `median_dolp_abs = 0.0922075` (claimed
0.0055043, gate threshold 0.02 — FAILS), `p95_dolp_abs = 0.231859` (claimed 0.0286531 — FAILS).
Since the claim's own code state cannot produce the claimed values and the evaluation is
deterministic, the claimed IPRT numbers were never produced by the committed code+config+reference
combination. (STATUS_README/README quote them; the audit's M-2 "unverifiable claims" risk is now a
demonstrated contradiction for this family.)

Physics localization (all probes in `results/validation/benchmark_iprt_a1_vector*_comparison.csv`
and `validation_report__split_iprt*.json`):
1. **Single scatter is exact.** `single_scatter_only=true` variant (`benchmark_iprt_a1_vector_ssonly.cfg`):
   model DoLP equals analytic Rayleigh sin^2/(1+cos^2) to 4 decimals at every vza (0.9415 vs 0.9415
   at vza 80), and model absolute intensity matches the analytic single-scatter slab integral to
   0.1–0.8 % in F0=1 units.
2. **The defect is the polarization of the multiple-scatter components.** Decomposition at vza 80:
   reference MS (MYSTIC total minus analytic ss) = 0.0344 with implied MS DoLP ~ **+0.48 (aligned)**;
   model MS = 0.0296 (-14 % in intensity, acceptable) but implied MS DoLP ~ **-0.22 (ANTI-aligned)**.
   A uniform ~x0.6 DoLP deficit at all angles follows. This is a sign/frame error inside the
   deterministic-second and/or MC higher-order polarized accumulation — same bug family as the
   Marseille output-frame reflection, but at a DIFFERENT layer (the output-stage `incoming_light`
   transform cannot touch it; total Q sign here is still correct).
3. **Signed diagnostics** (build_fix, first production use): `diag_q_sign_agreement_frac = 1.0`
   (no output-frame flip in the benchmark meridian basis), `diag_signed_un_median_abs = 0.0016`
   (U essentially perfect), `diag_signed_qn_median_abs = 0.0922` (deficit is pure Q magnitude).
4. **Coherence with Marseille:** twilight Marseille is ss-dominated with MoM suppressing
   higher-order (raising DoLP -> observed +0.134 over-polarization); IPRT daytime tau=0.5 has
   anti-polarized MS (lowering DoLP). Both point at the second/higher-order polarized machinery.
   Next discriminator (Phase B): separate deterministic-second vs MC-higher contributions to the
   anti-alignment (e.g. a variant disabling deterministic second scatter, or capping MC orders).

## Operational notes added 2026-07-29

- compute1 MaxTime = 3 d, **compute2 MaxTime = 10 d** (use compute2 for anything > 3 d).
- The 1x1-bin wrapper prefix is effectively SINGLE-THREADED (OpenMP is over directions) and the
  46-wavelength twilight prefix at 384 photons did not finish a 10-min login probe; narrowing the
  wrapper top-level to one wavelength made it seconds. Wavelength range in the wrapper does not
  leak into cases (fresh load, and both benchmark cfgs pin 550 nm themselves).
- The frozen `_batched_work` shard checkpoints all carry mtime 2026-04-26 (bulk copy) — historical
  per-direction timing is NOT recoverable from mtimes.
- `git worktree add <scratch> 232f67f` + `-DCMAKE_CXX_STANDARD_LIBRARIES=-lstdc++fs` (GCC 8) is the
  recipe for claim-era rebuilds; worktree removed after use.

## New artifacts 2026-07-29

| Path | What |
|---|---|
| `config/validation_split_{disort,iprt,rozenberg,koomen}.cfg` | single-case gate wrappers |
| `config/benchmark_iprt_a1_vector_ssonly.cfg` + `config/validation_split_iprt_ssonly.cfg` | single-scatter IPRT diagnostic variant |
| `config/paper_cases/..._measurement_branchcap_subset_fixed.cfg` | Phase B5 confirmation config (frozen budget + fixes) |
| `slurm/phase_a_split_cases.slurm`, `slurm/phase_a_split_koomen.slurm`, `slurm/validation_default_gate.slurm`, `slurm/branchcap_subset_fixed.slurm` | replacement jobs 779980–779983 |
| `results/validation/validation_report__split_disort{,_probe_login}.json`, `__split_iprt.json`, `__split_iprt_diag_buildfix.json` | landed artifacts incl. signed-QU diagnostics |
| `results/validation/benchmark_iprt_a1_vector{,_ssonly}_comparison.csv` | per-point IPRT model-vs-MYSTIC tables behind the finding |
| `tools/reconcile_validation_numbers.py` (extended) | label-aware family reconciliation (convergence restricted to the full gate) |

## Runtime confirmation of the explosion + gate-case exposure closed (2026-07-29, second concurrent session)

Written by a second session that diagnosed the 778296/778297 timeouts independently from the
runtime side; append-only, no edits to sections above. Cross-check against the config-side
analysis above — the two arrive at the same root cause.

**Live-process proof.** A resumed `strict_subset_fixed` batch_0000 worker (login node and a
compute1 probe on c065, job 779978, since scancelled) runs at 100 % CPU, state R, with zero
checkpoint advance; gdb stacks (3 snapshots, 20 s apart) show a `traceBandPath` self-recursion
tower (two recursive call sites, 0x431a71/0x431d9f in `build_fix`) grinding
`opticalDepthToSun -> computeOpticalProperties` exp/pow at every level. Not blocked — computing
a combinatorially exploding walk. Sidecar on c065 (`/work/vmo703/scratch/ls_probe_778297/
compute_probe-779978.out`): marseille worker ~100 % of a core, gcc8 ValidationRunner ~2 cores,
MemAvailable flat. Environment theories eliminated: sacct shows ZERO co-tenant jobs on c030 and
c035 for the entire 48 h window; nodes healthy, no reboot.

**778296 forensics.** Zawada single+multiple artifacts LANDED and pass every physics metric
(their `exit=1` is only the `measurement_case_present` quirk — smoke cfgs set no
`measurement_case_config`). The default gate then entered `default_clear_sky.cfg`, whose
`measurement_case_config=measurement_rozenberg_hminus6.cfg;measurement_koomen_meridian_
hminus6_polarization.cfg` reaches the same SZA=96 explosion (convergence full-sky may have
legitimately used the first hours; the twilight measurement cases are non-terminating).

**Claim-integrity finding extended to Koomen/Rozenberg.** Both case cfgs enter git at 232f67f
(2026-03-31), the same commit that introduces `twilight_higher_order_guiding` (default TRUE,
`twilight_higher_order_branches` default 2) with NO cap knob anywhere —
`higher_order_recursive_branch_cap` first appears at f05e47d (2026-04-30). The 232f67f
`traceBandPath` already has the same recursive branch loop. Therefore the SZA=96 all-orders
Koomen/Rozenberg cases are NON-TERMINATING at their own claim-era commit, and any historically
claimed numbers for them cannot have come from completed runs of the committed code+config —
same contradiction family as the IPRT finding above, established from the runtime side.

**Actions taken (2026-07-29 ~19:30Z).**
1. Appended `higher_order_recursive_branch_cap=1` (+ provenance comment) to
   `config/measurement_koomen_meridian_hminus6_polarization.cfg` and
   `config/measurement_rozenberg_hminus6.cfg`. Declared deviation: claim-era defaults are
   non-terminating; cap=1 is the marseille production-tier semantics and the only runnable
   setting. Configs are read at case start, so the still-PENDING jobs 779981 (koomen),
   779982 (full gate), 779983 (marseille branchcap) inherit the fix without resubmission.
2. Job 779980 had already entered `split_rozenberg` at 18:53:34Z with the uncapped cfg ->
   scancelled at ~22 min elapsed and resubmitted as **780007** (disort/iprt re-land in ~17 s
   each; rozenberg now capped).
3. Retired diagnostics: freeze-probe job 779978 scancelled; login probes killed. Probe evidence
   preserved under `/work/vmo703/scratch/ls_probe_778297/`.
4. Stale artifact note: `_batched_work/*strict_subset_fixed*` (16 batch cfgs + band-1
   checkpoints) were materialized from the uncapped tier and are superseded by the branchcap
   tier (779983, distinct case_id, no filename collision). Clean up at owner's discretion.

**Job ledger after this session's actions:** 779981 koomen (PENDING, capped via patch),
779982 full gate (PENDING, twilight cases capped via patch), 779983 marseille branchcap
(PENDING, cap in cfg), 780007 split disort/iprt/rozenberg (PENDING, rozenberg capped via patch).

## 780007 outcome + Rozenberg reconciliation (2026-07-29 ~19:45Z)

**780007 completed all three stages in 52 s** (disort 17 s, iprt 17 s, rozenberg 18 s with
`higher_order_recursive_branch_cap=1` — from non-terminating to 18 seconds; sacct shows FAILED
only because the driver propagates the cosmetic `overall_pass=false` exit code; all artifacts
archived). disort/iprt re-landed bit-identical to the 779980/login values.

**Rozenberg claim contradicted:** fresh `normalized_rmse = 0.264325` vs claimed `0.0849321`
(ratio 3.1, FAILS the 0.1-class gate expectation implied by the claim). Framing matters and
differs from the IPRT case: the claim-era config was NON-TERMINATING (finding above), so the
claimed number cannot have come from any completed run of the committed config; the fresh number
is the first real number for this case and embeds the declared cap=1 deviation. Koomen (779981,
same family, same expectation of contradiction) will complete the picture.

Score after this event: claims REPRODUCED for DISORT (5 digits) + Zawada x2 (MC noise);
claims CONTRADICTED for IPRT (16.8x on DoLP) + Rozenberg (3.1x on RMSE); pending: Koomen,
convergence. STATUS_README/MODEL_ASSUMPTIONS/README corrections (Phase A5) should be drafted
once Koomen lands, quoting fresh numbers with the cap deviation declared.

## Koomen reconciliation + default gate started (2026-07-29 ~20:45Z)

**779981 (Koomen split, capped) completed in 19 s.** Verdict MIXED — unlike Rozenberg/IPRT:
`normalized_rmse` fresh **0.102507** vs claimed 0.0470789 (2.2x — the claim remains unproducible
as configured, same non-terminating-family argument as Rozenberg, but the fresh value still
PASSES the gate threshold), while `p95_dolp_abs` fresh **0.0309199** is BETTER than the claimed
0.0366844. So the Koomen case is scientifically healthy under the declared cap deviation even
though the historically quoted RMSE was never real.

**779982 (full default gate) started 20:42:33Z** on compute2 — the multi-day convergence stage
is now the only outstanding Phase A artifact besides the branchcap subset (779983, still queued).

Final Phase A claim scoreboard (all four cases + Zawada now measured):
| Family | Verdict | Fresh vs claimed |
|---|---|---|
| DISORT scalar | reproduced, 5 sig figs | 1.000 |
| Zawada single + multiple | reproduced, MC noise | 0.99-1.19 |
| Koomen | mixed: RMSE 2.2x claim but passes; DoLP better than claim | 2.18 / 0.84 |
| Rozenberg | contradicted, FAILS | 3.1x |
| IPRT A1 vector | contradicted, FAILS (real MS-polarization defect) | 3.9x / 16.8x |
| convergence | pending 779982 (started) | — |

## 779983 branchcap-subset confirmation LANDED — Phase B5 closed (2026-07-30, second session)

Correction to the 20:45Z entry above: 779983 was not "still queued" — it COMPLETED 2026-07-30
01:45Z (started 01:42Z), two minutes before that entry was written. 48/48 directions, 256/256
samples each, **209 s total runtime on 32 workers** (budgeted 6 days). That is the cap=1 collapse
demonstrated in production: the identical workload uncapped did 0/48 directions in 48 h.

**Scoring vs the subset-specific predictions (all three hit):**

| Metric (like-for-like) | strict_subset era | Predicted | Actual 779983 | Verdict |
|---|---|---|---|---|
| raw median \|AoP err\|, all 48 rows | 52.5 deg | 21.2 deg | **22.29 deg** | HIT (within 1.1 deg) — frame fix confirmed on subset |
| masked 31/48 peak-norm shape RMSE | 0.2531 (MoM, r2-derived) | 0.2554 | **0.2085** | better than predicted |
| signed DoLP err median | — | +0.134 persists | **+0.1431** (mean +0.1303) | persists — genuine physics gap stands |

Metric-definition notes (needed to avoid false mismatches when quoting the report file):
- The report's `normalized_rmse=0.20852` IS the masked metric: mask = `normalized_reference >= 0.05`
  keeps exactly 31/48 rows and reproduces the value to 4 decimals (verified from the comparison CSV).
- The report's `median_aop_deg=32.23` is the DoLP-masked AoP variant (reference DoLP above threshold
  only); the prediction's "raw median" over all 48 rows is 22.29. Quote whichever matches the metric
  being compared, not interchangeably.
- Shape RMSE beating the r2-derived prediction (0.2085 vs 0.2554) is not decomposed: candidates are
  fresh RNG realization and the live fix-stack/quadrature knobs vs the r2-field re-aggregation the
  prediction was computed from. Not load-bearing for any pass/fail call.

Region detail for the DoLP-gap physics work: solar-vertical midzen mean signed bias +0.307
(second_frac 0.56, second_rr_frac 0.67 dominant), antisolar midzen +0.105, bright-horizon arc
+0.085, flip-215 sector is the only negative region (-0.137). Over-polarization concentrates
where deterministic second-order RR dominates — consistent with the IPRT MS-polarization defect
being the same family.

**Archived to `paper/tables/reproducibility/`**: the 779983 trio (report .txt, comparison.csv,
region_summary.csv), all four split validation reports + their comparison CSVs (780007/779981),
and `reconcile_snapshot_2026-07-30.txt` (label-aware reconcile output, 6/7 families fresh).

**Outstanding**: only 779982 (full default gate) — started 01:42:33Z 07-30 on c103, ~18 h into the
convergence stage as of this entry, 9 d walltime, twilight measurement cases capped via the cfg
patches. Phase A5 claim-correction drafting (STATUS_README/MODEL_ASSUMPTIONS/README) is unblocked
now that Koomen has landed; only the convergence family still awaits fresh numbers.

## Phase A5 claim corrections APPLIED to the three quoting docs (2026-08-01)

The reconciled numbers are now in the docs themselves, not just in this log:

- `STATUS_README.md` "Validation status": the "everything passes" framing paragraph replaced with
  a pointer to the revalidation; the benchmark-gate and measurement-gate bullets corrected; the
  "Current passing validation metrics" block replaced by a per-family reconciled block (fresh
  value, threshold, pass/fail, historical quote) with the IPRT unreproducibility + MS
  anti-polarization finding and the cap=1 deviation declared; the provenance sentence rewritten
  (split reruns for the four cases, April-2026 Windows values only for convergence, pending
  779982); one new bullet in "What is still blocking paper use" tying the failures + the
  Marseille +0.14 DoLP bias into the blocking list. Zawada smoke bullets annotated "(reproduced
  within Monte Carlo noise in the 2026-07-27 Linux rerun)".
- `monte_carlo_cpp/README.md` "Current validation result": convergence bullet re-scoped to the
  April 2026 Windows run pending 779982; scalar/vector benchmark bullet and measurement bullet
  corrected with fresh values + thresholds + the cap deviation; Zawada bullets annotated.
- `MODEL_ASSUMPTIONS.md`: "default passing benchmark/measurement suite" wording de-passified and
  a dated revalidation-status bullet added to each; the "still missing before paper-safe use"
  bullet about vector polarization rewritten — vector-polarization support currently rests on the
  Zawada smoke benchmarks + analytic single-scatter checks, and the default gate presently has no
  passing multiple-scatter vector benchmark.

**New fact surfaced while pulling thresholds** (updates the 07-29 scoreboard's "Koomen MIXED but
passing" wording): fresh `koomen median_dolp_abs = 0.0309199` **fails** its `0.02` threshold in
`validation_report__split_koomen.json`. That metric is not among the historically quoted values
(only RMSE and p95 were quoted), but under the current gate definition the Koomen case fails its
median-DoLP criterion — so Koomen is "mixed: RMSE passes at 2.2x the historical quote, p95 better
than quoted, median-DoLP fails", not a clean pass. All three docs state it that way.

Convergence rows in all three docs are explicitly marked "historical, pending job 779982"; when
that report lands (`validation_report__default_clear_sky.json`), update those rows + this log and
archive the report to `paper/tables/reproducibility/`. Job watcher re-armed 2026-08-01 (monitor
bs8udi9lh: state transitions, gate-stdout growth, report landing; 300 s poll) after the previous
monitor died with the session.

## FOUND AND FIXED: the multiple-scatter anti-polarization is a chi-sign handedness bug in eventMuellerMatrix (2026-08-01)

Task #19 executed end-to-end in this session. Chain of evidence:

**1. Order decomposition localized it to deterministic second order.** A standalone
`MonteCarloCPP config/benchmark_iprt_a1_vector.cfg` run dumps per-bin order-resolved Stokes
(`results/benchmark_iprt_a1_vector/sky_result.csv`). With the frames sign-matched via the exact
single scatter (model `first_Q` and MYSTIC `q` are both negative), the MS residual splits as:
MYSTIC MS q~ runs -0.03 -> -0.44 with VZA (single-scatter orientation); the model's deterministic
second order (pure RR here) runs **+0.03 -> +0.46 — sign-flipped with |q~| tracking the
single-scatter magnitude** (45 deg: +0.319 vs SS -0.333); the MC higher orders were statistically
zero (|Q|/sigma ~ 0.04 at 16384 photons). Analysis script:
`scratch session file iprt_ms_decomp.py` (rerunnable; takes the sky CSV path as argv[1]).

**2. A dipole-ground-truth probe proved the mechanism.** `det2_frame_probe.cpp` (session
scratchpad) replicates the exact producer->consumer chain
(`deterministicSingleScatterAlongRayBreakdownAllBands` -> det2 `eventMuellerMatrix` application)
and compares against an independent E-vector/dipole double-scatter computation over 45 geometries:
coplanar cases match exactly; every out-of-plane case is wrong in Q and U (plus up to 0.2% in I —
an output-side rotation error could not touch I). Cause: both chi angles are measured about the
BACKWARD ray axes, but the chain applies `R(-chiIn)` (forward-propagation-convention rotation)
with `R(+chiOut)` (backward-convention rotation) — mutually inconsistent handedness. Unpolarized
inputs are insensitive to chiIn, so single scatter — the only polarized path the validation suite
constrains tightly — was correct all along, hiding the bug. Flipping EITHER rotation restores
consistency: variants V1 (`R(+chiIn)`) and V2 (`R(-chiOut)`) both match the dipole ground truth to
1e-9 in Q and I across all 45 geometries (they differ only in U handedness). V1 is the correct
production choice because it preserves the chiOut path that every validated output convention
(single scatter, Marseille AoP) is built on.

**3. Fix implemented, knob-gated, legacy bit-identical.** `event_frame_chi_sign_fix` (bool,
default false) added to the config struct, parser, and metadata JSON; `eventMuellerMatrix` takes a
`chiInSignFix` parameter applied at all 10 call sites (SS view path, SS-along-ray, det2 rr/ar/ra/aa,
MC NEE source event, MC scatter matrix). With the knob off the new binary reproduces the previous
IPRT sky CSV **bit-identically**, and all 5 ctest suites pass — frozen-artifact reproducibility is
untouched. Rebuilt in `build_fix`.

**4. With the fix on, the "unproducible" IPRT claims reproduce almost digit-for-digit.**
`validation_split_iprt_chifix.cfg` -> `benchmark_iprt_a1_vector_chifix.cfg` (RUNNER must point at
`build_fix/ValidationRunner`; the wrapper script defaults to `build_current`, which silently
ignores the knob — first attempt produced bit-identical legacy values because of exactly that):

| metric | chifix fresh | claimed 2026-04 | legacy fresh | gate |
|---|---|---|---|---|
| median_dolp_abs | **0.0055043** | 0.0055043 | 0.0922075 | PASS (thr 0.02) |
| p95_dolp_abs | 0.0286162 | 0.0286531 | 0.231859 | PASS (thr 0.05) |
| median_intensity_rel | 0.00546554 | 0.00561889 | 0.0220254 | PASS |
| p95_intensity_rel | 0.0178296 | 0.0176791 | 0.0710573 | PASS |

Order decomposition with the fix: model MS q~ median -0.2049 vs MYSTIC -0.2104; det2 flips to
-0.32; MC higher orders now carry real negative polarization (-0.11 median). The residual
high-VZA gap (-0.33 vs -0.44 at 75 deg) is bin-center-vs-point smearing + MC noise, and the exact
benchmark evaluation above shows the actual pointwise metrics are claim-level.

**Reinterpretation of the Phase A finding.** "IPRT claims unproducible at their own commit"
stands as a statement about the RECORDED code (the from-source rebuild at 232f67f reproduces the
broken values bit-identically). But the claims themselves are now vindicated: they match the
consistent-chain physics to 4+ significant figures, so the claim-era numbers were evidently
produced by a local/unrecorded code state that had the consistent rotation, and the committed
history carried the sign regression from the start. The same reinterpretation may apply to
Rozenberg/Koomen (their claims are also from that era); chifix reruns of both are in flight.

**In flight at the time of writing:** `validation_split_zawada_chifix` (does the currently-passing
Zawada pair survive/improve?), `validation_split_rozenberg_chifix` + `validation_split_koomen_chifix`
(do the failing/mixed twilight cases collapse toward their claims?), and Slurm job **781658**
(`branchcap_subset_chifix.slurm`): the 48-direction Marseille subset with the fix on — baseline
779983 signed DoLP err median +0.1431, solar-vertical midzen bias +0.307 concentrated where
second_rr dominates; if the chi bug is the Marseille over-polarization mechanism, that number
drops. Production/frozen configs remain on legacy (knob off) pending a deliberate adoption
decision (Phase E).

## Chi-fix blast radius complete (2026-08-01, later): Zawada confirms it, twilight gap is a separate mechanism

All follow-up runs landed the same day. Full verdict table for `event_frame_chi_sign_fix=true`
(all runs via `run_phase_a_validation.sh` with `RUNNER=build_fix/ValidationRunner`, reports
`validation_report__split_*_chifix.json`; Marseille via Slurm job 781658, COMPLETED 48/48 in
3 m 31 s):

| case | legacy (fix off) | chifix (fix on) | reading |
|---|---|---|---|
| IPRT A1 vector | FAILS DoLP (0.0922/0.2319) | **PASSES all four; equals claims to 4+ sig figs** | claims vindicated; committed code had the regression |
| Zawada single (SS-only) | passes | **bit-identical** | knob is a proven no-op on unpolarized-input paths |
| Zawada multiple (all orders) | passes (DoLP 0.0199/0.0354) | **passes, 2.8-3.6x better (0.00706/0.00981)**; intensity also improves | independent published benchmark prefers the fixed chain |
| Rozenberg twilight intensity | fails 0.264325 | fails 0.26407 | unaffected (intensity-dominated); claim stays unexplained |
| Koomen twilight polarization | RMSE 0.1025 pass; DoLP 0.0309 | RMSE 0.0999 pass; **DoLP worsens to 0.0515 (p95 now fails too)** | twilight over-polarization deepens under correct rotations |
| Marseille 48-dir subset | signed DoLP med +0.1431 (mean +0.1303); raw AoP med 22.29; RMSE 0.2085 | **+0.1938 (mean +0.1745)**; AoP 23.37; RMSE 0.2070 | DoLP gap is NOT the chi bug — it widens; solar-vertical midzen bias +0.307 -> +0.401 |

Operational trap logged on the way: the first Zawada chifix attempt used the
`benchmark_zawada_spherical_vector_*_smoke.cfg` files and "exploded" (median_intensity_rel 0.52,
p95 = 1.0). A knob-off control through the identical wrapper produced bit-identical garbage: the
`*_smoke.cfg` files are a stale variant tier that does not match the bundled reference (the
validated Jul-27 runs and the STATUS_README "smoke" quotes actually evaluate the BASE
`benchmark_zawada_spherical_vector_{single,multiple}.cfg` — their metric prefixes carry no
"smoke"). The base-case legacy control through my wrapper reproduces the archived Jul-27 numbers
at ratio 1.0000 on all eight metrics. Do not use the `*_smoke.cfg` case files for anything.

Standing interpretation after all of this:
- The chi-sign fix is correct (dipole ground truth, IPRT claim match, Zawada-multiple external
  improvement, exact no-op where it must be one) and should be the production default for any
  non-frozen future work; frozen tiers keep knob-off reproducibility.
- The Marseille/Koomen twilight over-polarization is a genuinely separate physics gap (finding 3),
  now sharpened: it persists and worsens under a correctly-rotating chain, so its mechanism is
  missing depolarization physics (candidate territory: aerosol/molecular depolarization terms,
  ocean/ground and horizon source terms, spectral band coverage), not frame handling.
- Rozenberg (intensity) and the Rozenberg/Koomen claim-era numbers remain unexplained and
  unproducible; the cap=1 deviation framing in the docs stands.

Chifix artifacts archived to `paper/tables/reproducibility/`: the five split chifix reports (+
zawada base legacy control), the IPRT chifix comparison CSV, and the Marseille chifix trio
(report .txt, comparison csv, region summary).

## 779982 diagnosed as recursion-doomed; capped-deviation gate submitted as 783441 (2026-08-03)

At 4 d 19 h of its 9 d wall, 779982's `default_clear_sky_stdout.log` is still the 523-byte
startup burst — zero stage completions. Root-cause check: `config/default_clear_sky.cfg` carries
NO `higher_order_recursive_branch_cap` line, so the gate's own convergence skies (19x36 = 684
directions x 46 bands, all orders) run the SZA-97 twilight sun under the hpp defaults
(guiding on, branches 2, **cap 0**) — the same runtime-confirmed non-terminating regime that ate
778296 and the strict subset. Cost scaling from 779983 (~140 s/direction single-thread at cap=1,
25 bands, 256 photons) says a CAPPED version of both convergence skies is ~2 h at 36 threads —
so 4.8 d of silence is diagnostic, not slow-but-fine. The 9 d budget in
`validation_default_gate.slurm` was extrapolated from 778296's 48 h as if that had been real
progress; it was recursion spin, so the extrapolation was circular. Note this also means the
historical convergence claims (0.0382187/...) are in the same evidentiary bucket as the IPRT
claims: not producible by the recorded config on this machine class (the sky would never finish),
so they too presumably came from a claim-era local state (pre-guiding code, or different local
knobs, or the unrecorded consistent-chain build).

Action (2026-08-03): left 779982 running to walltime as the legacy control (canceling a 4.8 d
investment on inference alone is not warranted; if it defies the analysis and completes, its
report is strictly better evidence), and submitted **783441** `validation_default_capped.slurm`
(compute2, 36 cpu, 2 d wall): the full gate with `config/validation_default_capped.cfg` = the
default cfg + `higher_order_recursive_branch_cap=1` as a DECLARED deviation (in-file provenance
comment; same deviation class as the Rozenberg/Koomen patch, and it matches the frozen production
tier, so its `convergence_*` metrics are the production-configuration self-convergence numbers).
Report will land as `validation_report__default_clear_sky_capped.json`; reconcile the convergence
family against the historical quotes from there, with the deviation declared, and update the
convergence rows in STATUS_README / monte_carlo_cpp README / MODEL_ASSUMPTIONS (currently marked
"pending 779982" — they should end up saying "uncapped claim configuration nonterminating;
capped-deviation revalidation gives ..."). Monitor bvabu0a6x watches both jobs (state changes,
stdout growth, either report landing; the old shared `validation_report.json` mtime watch was
dropped because every local wrapper run false-fired it).

## 783441 landed: convergence family PASSES under the declared cap=1 deviation — reconciliation complete, tasks #12/#13 closed (2026-08-03)

Job 783441 finished the full default gate **in 21 seconds** (Slurm state FAILED = expected
nonzero exit from `overall_pass=false`; all four case stages + convergence ran and the report
landed as `validation_report__default_clear_sky_capped.json`). The 21 s wall is itself the
sharpest possible confirmation of the 779982 diagnosis: same Apr-27 `build_cluster_gcc8`
binary (has the cap knob — verified via `strings`), same config except
`higher_order_recursive_branch_cap=1`, and the whole gate costs 21 s — while the uncapped
779982 sat 4 d 19 h without completing a single stage. My earlier "~2 h" estimate extrapolated
from the Marseille per-direction cost was off by ~250x (the default sky's pure-Rayleigh
atmosphere is far cheaper per direction than the aerosol-laden Marseille tier); the estimate
being wrong makes the diagnosis STRONGER, not weaker. Verified SZA from the recorded config:
cos(SZA) = cos(45°)·cos(100°) → SZA ≈ 97.05°, so the recorded gate's convergence skies are
squarely in the runtime-confirmed non-terminating uncapped-twilight regime.

Fresh convergence family (coarse 128 vs fine 256 photons/bin, 19x36 dirs x 46 bands,
`evaluateConvergence` full-sky pairs):

| metric | fresh (capped) | threshold | historical claim | ratio |
|---|---|---|---|---|
| convergence_peak_intensity_rel | 0.0251017 | 0.05 PASS | 0.0382187 | 0.66 |
| convergence_peak_dolp_abs | 0.0146482 | 0.02 PASS | 0.00138362 | 10.6 |
| convergence_flux_rel | 0.00645125 | 0.02 PASS | 0.00211343 | 3.05 |

All three pass; the values differ from the April 2026 quotes as expected — the recorded
uncapped config can never terminate, so the historical numbers are from an unrecorded
claim-era local code state (same bucket as the IPRT claims). Internal consistency check: the
capped report's DISORT/IPRT/Rozenberg/Koomen stage values reproduce the 2026-07-29 split-rerun
values digit-for-digit (0.0188114 / 0.0922075-legacy / 0.264325 / 0.102507+0.0309199),
confirming identical code paths and the case_id-independent RNG seeding.

Reconciliation tooling: `tools/reconcile_validation_numbers.py` gained
`default_clear_sky_capped` as a convergence-family fallback label and
`split_zawada_base_legacy` as the preferred zawada label (the `*_smoke` tier is stale). Full
run now reports **0/7 families missing** for the first time; output archived as
`paper/tables/reproducibility/reconcile_validation_numbers_2026-08-03.txt` alongside the capped
report JSON. Doc rows updated in all three docs (STATUS_README "Validation status" convergence
bullet + provenance line; monte_carlo_cpp/README validation bullet; MODEL_ASSUMPTIONS
Validation Gates bullet). Tasks #12 and #13 closed.

Final claim-family scoreboard (legacy chain, declared deviations where noted):
convergence PASS (cap=1 deviation) · DISORT reproduced exactly · Zawada single+multiple
reproduced ratio~1.000 · IPRT fails legacy / vindicated under `event_frame_chi_sign_fix` ·
Rozenberg FAILS (0.264 vs claim 0.085, unexplained) · Koomen mixed (RMSE passes at 2.2x claim;
median-DoLP 0.0309 fails 0.02). 779982 stays queued to walltime as the legacy control
(expect TIMEOUT ~08-08 with an empty report); a fresh monitor watches only its exit.

## 779982 legacy control closed: TIMEOUT at 9 days, zero stages, no report (recorded 2026-08-09)

Final sacct record: `779982 TIMEOUT 9-00:00:12`, killed by the scheduler 2026-08-07T20:42:45.
Its entire stdout is four lines: driver start 2026-07-30T01:42:33Z, `=== run_case
default_clear_sky ===`, then the slurmstepd time-limit cancellation nine days later. No
`validation_report__default_clear_sky.json` exists. The prediction from the 08-03 diagnosis
held exactly: the recorded uncapped `default_clear_sky.cfg` (SZA 97.05 deg, no
`higher_order_recursive_branch_cap` line) is non-terminating, on the very binary
(`build_cluster_gcc8/ValidationRunner`) that completed the identical capped gate (783441) in
21 s. The uncapped-vs-capped control pair is now complete and archived in the docs; the three
docs' convergence rows were updated from "~5 days silent" to the final 9-day TIMEOUT fact.
Nothing further is pending in Phase A.

## Four-track kickoff: depolarization gap, Rozenberg, paper, chi-fix adoption (2026-08-09)

User greenlit all four post-Phase-A tracks in parallel (tasks #20-#23). Ordering constraint
honored: the chi-fix production FREEZE waits for the depolarization verdict so production
freezes once; everything else runs concurrently.

**Cap semantics caught before wasting jobs**: `higher_order_recursive_branch_cap` clamps
branchCount into [1, cap] at every order>=3 event and the ambient branch factor is 2, so
cap=2 IS the divergent ~2^61 regime — there is no cap ladder, cap is binary (1=linear chain,
>=2=non-terminating). The order-content probe is a `max_events_guard` ladder instead; the
runtime floors the guard at 8 (`std::max(8, guard)`), lowered to 2 in the new build_fix2
(bit-identical for every recorded config, all of which use guard>=8... i.e. 64).

**King-factor depolarization knob implemented** (`rayleigh_depolarization_factor`, default
0.0): depolarized Rayleigh matrix after Hansen & Travis (delta=(1-rho)/(1+rho/2), isotropic
completion keeps f11 normalized, f44 carries delta'=(1-2rho)/(1-rho)). Enters ONLY via
eventMuellerMatrix (all 10 call sites) — samplers/pdfs stay pure as importance-sampling
proposals, so the estimator remains unbiased. build_fix2: ctest 5/5; analytic probe
DoLP(90)=(1-rho)/(1+rho) to 1e-16, f11 integral exactly 1, rho=0 bit-identical to the pure
matrix. Runtime bit-identity smoke (build_fix vs build_fix2 on the tiny_smoke case) running
on the login node; the build_fix2 Slurm leg (guard 3/4/6 + rho=0.0279 subset runs,
slurm/depol_lowguard_fix2.slurm) queues after it passes.

**Guard ladder result (job 787249, 12 min)**: on the 48-direction Marseille subset with the
corrected chain, signed DoLP median is +0.1906/+0.1938/+0.1938/+0.1938 at guard
8/16/32/64 — guard16=guard32=guard64 exactly, guard8 barely different, higher_frac median
0.062-0.072 throughout. The MC chains effectively die by order ~10-16 and the
over-polarization is fully formed at LOW orders: it is not a deep-order truncation artifact.

**Rozenberg profile finding (task #21, from existing split-run CSV, no new compute)**: the
model's normalized sky falls off the horizon peak 3-15x FASTER than the Rozenberg reference
at every zenith/azimuth (mod/ref 0.38 at zen80 down to ~0.07-0.10 over much of the sky), and
the model is single-scatter dominated (first_frac 0.4-0.95) where the real SZA-96 twilight
sky is multiple-scatter dominated. The median_of_means suspect is ELIMINATED (neither
Rozenberg nor Koomen cfg sets robust_groups, so the estimator is already plain mean). Working
hypothesis: the twilight MS field is underproduced in absolute terms at all depths — one
shared root cause for the intensity-profile collapse AND the +0.19 over-polarization
(too little MS = too little depolarized light). Discriminator submitted: job 787252
(`validation_split_rozenberg_noguide`, twilight_higher_order_guiding=false = vanilla
single-branch phase-sampled continuation, unbiased control): if its MS yield is materially
larger, the guided cap=1 estimator is biased low; if unchanged, the deficit is physics
(refraction? aerosol profile?).

**Chi-fix candidate gate — two findings (job 787250, 21 s)**:
1. Wrapper-knob trap: case stages load their cfgs FRESH, so the wrapper's
   event_frame_chi_sign_fix reached only the convergence stage (case stages bit-identical to
   legacy). Candidate gate re-pointed at the *_chifix case variants and resubmitted (787253).
2. NEW ANOMALY: under the corrected chain the default gate's convergence family FAILS —
   convergence_peak_intensity_rel 0.290216 (thr 0.05; legacy-capped 0.0251), flux_rel
   0.0239651 (thr 0.02), peak_dolp 0.0196 (marginal pass). The corrected chain appears to
   have much heavier-tailed higher-order intensity at the gate's plain-estimator
   256-photon settings. MUST be diagnosed before any production freeze (adoption blocker).

**Paper (task #22, done)**: `paper/DRAFT_OUTLINE_2026-08-09.md` — framing A (frozen
calibrated-pipeline story) vs framing B (model + validation + regression-archaeology,
recommended), full section skeleton, figure inventory F1-F6, claims discipline, blocking
inputs. Frozen April package untouched.

**In flight at close**: 787251 full-field 683-direction corrected-chain Marseille run
(adoption numbers + region-resolved depol data), 787252 Rozenberg no-guiding control,
787253 corrected candidate gate, login-node bit-identity smoke, then the build_fix2
guard3/4/6+depol0279 leg.

### Same-evening results: candidate gate v2 + Rozenberg no-guiding null (2026-08-09)

**787253 candidate gate v2** (case stages re-pointed at *_chifix cfgs, 22 s): IPRT chifix all
four PASS at the vindicated values (median_dolp_abs 0.0055043 exact); Koomen chifix RMSE
0.0998924 pass (slightly better than legacy 0.1025) but DoLP 0.0515364 fails both thresholds;
Rozenberg chifix 0.26407 fail (untouched); and the convergence anomaly REPRODUCES EXACTLY
(peak_int 0.290216 — deterministic, so it is a stable property of the corrected chain at the
gate's 256-photon plain-estimator settings, not a flake). Report:
validation_report__default_clear_sky_capped_chifix.json.

**787252 Rozenberg noguide control** (18 s): normalized RMSE 0.257042 vs guided 0.264325 —
the profile does not move without guiding. The guided cap=1 continuation is NOT suppressing
the multiple-scatter yield; second estimator suspect eliminated. The MS deficit is a
physics/inputs question: leading candidates are atmospheric refraction (Earth-shadow height
at SZA 96), the h^-6 aerosol profile realization, spectral weighting vs the Table-1
reference wavelength, and the TOA cap. (Wrapper oddity noted: this clone reports
measurement_polarization_reference_present=1/FAIL where the original split wrapper reported
0/PASS — informational metric, does not affect the RMSE readout; check the flag's semantics
if the wrapper is reused.) Report: validation_report__split_rozenberg_noguide.json.

### Full-field corrected chain landed + build_fix2 sealed (2026-08-09 evening)

**787251 full-field 683-direction corrected-chain Marseille run** (55 min, report
frozen_marseille_twilight_..._measurement_full_chifix.*): signed DoLP bias median **+0.2068**
(mean +0.1831, p95 +0.454), median |AoP err| 17.01 deg, higher_frac median 0.072. The
zenith structure is the fingerprint of the deficit hypothesis: signed bias median **+0.350
at zenith 0-30 deg, +0.142 at 30-60, +0.093 at 60-90** — the over-polarization is largest
exactly where multiple scattering should dilute DoLP the most (zenith), and smallest where
single scatter legitimately dominates (bright horizon). Consistent with the July B1 raw-field
finding and with the Rozenberg profile deficit: one under-produced twilight MS field explains
both. This run is also the chi-fix adoption candidate full-field artifact (task #23).

**Bit-identity trap + seal**: the first bitcheck ran the tiny_smoke config as-is — which is a
PRE-CAP uncapped SZA-96 config (zero cap lines) and therefore non-terminating (the 779982
trap in miniature; 111 CPU-min before I killed it). With cap=1 appended the case runs in
0.1 s, and **build_fix vs build_fix2 are bit-identical on all physics report lines and the
full comparison CSV** (only wall-clock timing_* telemetry differs). RULE reinforced: any
SZA>90 all-orders config MUST carry higher_order_recursive_branch_cap=1 — including scratch
copies for smoke tests. The fix2 leg (guard3/4/6 + depol0279) submitted as job 787262.

## Depolarization-gap decomposition COMPLETE: order depth is not the lever; King factor is ~14%; the rest is the missing MS bulk (2026-08-09 night)

Job 787262 (fix2 leg, 14 min) completed the experiment set. Full ladder on the 48-direction
corrected-chain subset (signed DoLP bias median / mean / higher_frac median):

| case | sgn med | sgn mean | hi_frac | sgn med (zen<45) |
|---|---|---|---|---|
| guard3  | +0.1925 | +0.1911 | 0.078 | +0.2515 |
| guard4  | +0.1784 | +0.1746 | 0.042 | +0.2519 |
| guard6  | +0.1925 | +0.1805 | 0.051 | +0.2724 |
| guard8  | +0.1906 | +0.1757 | 0.062 | +0.2619 |
| guard16 | +0.1938 | +0.1745 | 0.072 | +0.2617 |
| guard64 | +0.1938 | +0.1745 | 0.072 | +0.2617 |
| depol0279 | +0.1705 | +0.1498 | 0.065 | +0.2327 |

**Verdict 1 — chain depth is irrelevant**: truncating the MC chain all the way down to
order ~3 leaves the bias at +0.18-0.19 (low-guard medians bounce within MC noise of a small
component; hi_frac 0.04-0.08 throughout). The model's existing higher-order field is too
small everywhere for its depth to matter. Combined with the guard16=guard64 exactness:
the over-polarization is carried by the order-1+2-dominated field.

**Verdict 2 — molecular (King-factor) depolarization is real but minor**: rho=0.0279 moves
the signed bias -0.023 (median), -0.025 (mean), -0.029 (zenith<45) — about **13-15% of the
gap**, matching the textbook expectation (DoLP(90) scales by (1-rho)/(1+rho)=0.9457). Right
direction, right magnitude, worth adopting as production physics in the next frozen tier —
but not the main story.

**Remaining gap ~+0.15-0.17**: the missing multiple-scattering BULK — the same deficit the
Rozenberg intensity profile shows at 3-15x, now quantified from the polarization side. The
twilight MS field itself is underproduced (hi_frac ~7% where the real SZA-96 sky is
MS-dominated); with both estimator suspects eliminated (median_of_means inactive; noguide
null), the mechanism hunt is now physics/inputs: refraction (Earth-shadow height at SZA 96),
aerosol profile/optical-depth realization in twilight geometry, spectral weighting of the MS
field. That hunt continues under task #21 (Rozenberg), which is now formally the same
investigation. Task #20 closed.

## ROZENBERG SOLVED (input archaeology) + convergence anomaly characterized (2026-08-09 night, session 2)

**Rozenberg: the claim is explained and the case passes — the recorded config points at the
wrong aerosol profile.** The case is named hminus6 after Rozenberg's Table-1 h^-6 aerosol
scenario, but `measurement_rozenberg_hminus6.cfg` loads the generic
`clear_sky_midlatitude.csv`, whose aerosol extinction falls far SLOWER aloft than h^-6
(1.2e-7 vs 1.5e-10 /m at 30 km). Along SZA-96 tangent paths (~10^6 m at 10-30 km) that
excess high-altitude aerosol over-attenuates the twilight illumination and collapses the
diffuse sky. With a true h^-6 profile (aerosol ∝ h^-6 above 5 km, anchored at the recorded
profile's 5 km value; `clear_sky_midlatitude_diag_hminus6.csv`):

    normalized RMSE 0.264325 -> 0.079095  (PASSES thr 0.15; historical claim 0.0849321)
    solar-vertical model/ref ratios: 0.38/0.16/0.21/0.09/0.16 -> 1.11/0.70/0.99/0.65/1.27

Third archaeology win (after IPRT code state and the non-terminating convergence config):
the claim-era Rozenberg run evidently used an h^-6-type profile that was never recorded.
Sensitivity map from the discriminator set (all ~20 s runs, reports archived under
validation_report__split_rozenberg_{noguide,noo3,noaero,aerox4,refrbound,trueh6}.json):
ozone-off lifts diffuse ratios 1.5-2x (tangent-path absorption is first-order);
refraction bound (SZA-0.57 deg) lifts them 1.3-1.5x (real, secondary); aerosol x4
restructures the whole field; no-guiding is a null (estimator exonerated).
NEXT BRIDGE (open): the Marseille frozen atmosphere profile may be similarly
aerosol-rich aloft -> same over-attenuation suppressing its MS field -> the +0.15-0.17
DoLP remainder. Probe = subset run with an h^-6-ized VARIANT of the frozen profile
(variant file only; frozen inputs untouched).

**Convergence anomaly: characterized, not photon-fixable.** Bin-level dumps
(scratch conv_diag set, MonteCarloCPP standalone) show the 0.290216 gate value is
peak-bin migration: at 128 photons one monster higher-order path (total path I-weight
~4.3e-3, constant across budgets — same seed stream) inflates bin (zen 87.6, az 325)
above the physical arch bin (zen 82.9, az 235; second_I = 4.29e-5 = 87% deterministic).
Photon ladder of the arch bin higher-order term under the corrected chain:
hi_I = 1.37e-6 / 4.71e-6 / 2.85e-6 / 8.40e-6 at 128/256/512/1024 photons, with
hi_sd/hi_I GROWING 6.5 -> 12.5 -> 14.7 -> 16.7 — the classic near-infinite-variance
signature; pairwise peak_rel: 0.290 (128v256), 0.038 (256v512), 0.101 (512v1024).
Legacy has the same disease milder (ratio 8.65 at 256; passed 128v256 at 0.025 partly by
tail luck). No unphysical Stokes anywhere (|QUV|<=I in all 684 bins). Interpretation: the
plain-mean higher-order estimator is heavy-tailed under the twilight-guided weights
(consistent with the July finding: rel SE median 62% at 256 samples), and the corrected
chain amplifies the tail (coherent Q alignment raises some guided paths' I-transfer).
Remedy options for the candidate tier gate (decision open): (a) convergence metric on the
deterministic components + flux (second_I is 87% of the arch and rock-stable), (b) masked
field-RMSE convergence instead of single-bin max, (c) budget bump alone is REFUTED.
Side implication for the MS story: a heavy right tail means typical realizations
understate the MS field at these budgets even though the estimator is unbiased —
second-order contributor to the deficit; the h^-6 input effect is the dominant term.

### Bridge probe result: h^-6 profile REFUTED for the Marseille DoLP gap — mechanisms split (2026-08-09 night)

Job 787269 (h^-6-ized frozen-profile variant, corrected-chain subset): signed DoLP bias
median +0.1938 -> **+0.2278** (WORSE; zenith<45 +0.2617 -> +0.2919; hi_frac 0.072 -> 0.055).
Consistent with Rozenberg once decomposed: stripping high-altitude aerosol brightens the
twilight INTENSITY field (attenuation effect — what fixed Rozenberg's normalized shape) but
removes depolarizing aerosol scattering events, leaving a purer-Rayleigh, MORE polarized
sky. Conclusion: **Rozenberg (intensity, input archaeology) and the Marseille DoLP bias are
different mechanisms.** The DoLP levers now stand: King factor -0.023 (real, adopt);
aerosol-profile substitution counterproductive; PROMOTED LEAD = sampling-realization
understatement of the MS field (conv_diag showed single paths carrying bin-mean-scale
weight; the July analysis put plain-mean rel SE at median 62% — typical realizations miss
the depolarized heavy tail). Direct test submitted: job 787270, subset at 8x budget
(photons 2048); if hi_frac climbs and the bias falls with budget, the remaining deficit is
sampling, and the production remedy is estimator/proposal redesign, not new scattering
physics.

### p2048 verdict + the f22 smoking gun (2026-08-09, closing the night)

**Sampling-realization: real but small.** Job 787270 (8x budget, photons 2048): median
higher_frac DOUBLES 0.0723 -> 0.1450 (the tail understatement is real — typical
realizations do miss MS content), but the signed DoLP bias moves only -0.009 median
(+0.1938 -> +0.1844; improves in just 28/48 directions). Extrapolated to convergence this
buys a few hundredths at most. Demoted to a secondary term.

**The structural finding: the model's aerosols cannot depolarize.** Both
`aerosol_phase_matrix.csv` files (Marseille frozen + optics reference) satisfy f22 = f11
EXACTLY at every wavelength and angle — the perfect-sphere Mie identity (f34=0 throughout
as well). An aerosol scattering event therefore transmits linear polarization losslessly;
the model's only depolarization channels are unpolarized surface reflection and (since
tonight) the molecular King factor. Real urban/marine aerosol has f22/f11 ~ 0.5-0.8 in
side/backscatter (nonsphericity), which depolarizes at every aerosol event — and aerosol
events are a large fraction of the twilight second order (second_ar/ra/aa fractions 0.24+
in the arch zones). This is the best-motivated remaining candidate for the ~+0.15 residual.

**Next experiment (queued for next session): `aerosol_depolarization_f22_ratio` knob**
(default 1.0 = legacy bit-identical; scales f22 relative to f11 at evaluation), subset scan
at 0.9/0.8/0.7/0.6 stacked with the King factor. Elimination table now: estimator choice x2,
chain depth, aerosol profile shape, budget realization — all eliminated/quantified small;
King -0.023 (adopt); aerosol nonsphericity = open prime suspect.

## Aerosol f22-ratio knob implemented + scanned: real but small; loading probe launched (2026-08-10)

Knob `aerosol_depolarization_f22_ratio` (default 1.0, bit-identical) scales f22/f33/f44 at
aerosol events inside eventMuellerMatrix, leaving f11 (energy/normalization) and f12
(polarization generation, hence all single scatter) untouched. build_fix3: ctest 5/5,
bit-identical to build_fix on the capped tiny case (build_fix2 left frozen as the 08-09 run
binary). Scan on the corrected-chain subset (job 787296, 18 min):

| case | sgn med | sgn(z<45) | AoP med | shape RMSE |
|---|---|---|---|---|
| baseline   | +0.1938 | +0.2617 | 23.37 | 0.1731 |
| king only  | +0.1705 | +0.2327 | 23.33 | 0.1738 |
| f22 0.9    | +0.1937 | +0.2630 | 22.60 | 0.1731 |
| f22 0.8    | +0.1910 | +0.2585 | 21.63 | 0.1731 |
| f22 0.7    | +0.1884 | +0.2553 | 20.70 | 0.1731 |
| f22 0.6    | +0.1857 | +0.2538 | 20.01 | 0.1731 |
| 0.7 + king | +0.1671 | +0.2246 | 20.63 | 0.1740 |

**Verdict: aerosol nonsphericity is real but SMALL for DoLP** (-0.0027 per 0.1 of ratio;
-0.008 at the aggressive 0.6) because the polarized twilight signal is dominated by
Rayleigh-Rayleigh chains the knob does not touch (second_rr_frac 0.56-0.70 in the arch
zones). Side findings: AoP median improves steadily 23.4 -> 20.0 deg (aerosol
depolarization cleans angle structure), and shape RMSE is bit-stable at 0.1731 across the
scan — the f11/f12-untouched design isolation confirmed empirically. Best combined physics
so far: 0.7+king = +0.1671 vs +0.1938 baseline (~0.027 of the ~0.19 recovered).

**Also notable**: the Marseille bias (+0.19) is 4-6x Koomen's DoLP error (0.031-0.052)
under the identical model chain — the two polarization references disagree about how wrong
the model is, which keeps a reference/reduction-side component on the Marseille table.

**Next probe in flight (787301)**: total aerosol loading x4 on the frozen Marseille profile
(variant file) — the one input never scanned on Marseille; a hazy Mediterranean August day
could far exceed the climatological profile, and even Mie aerosol depolarizes the SKY DoLP
by adding weakly polarized forward-scattered light.

## MECHANISM HUNT CONCLUDED: the Marseille DoLP gap is dominated by aerosol LOADING (2026-08-10)

Completing scan (787304: x2, x3, x3+King+f22r0.7) joins 787301 (x4). Full curve on the
corrected-chain subset (signed bias median / zen<45 / AoP med / shape RMSE / |bias| med):

| case | sgn med | z<45 | AoP | shape | |bias| |
|---|---|---|---|---|---|
| baseline    | +0.1938 | +0.2617 | 23.37 | 0.1731 | 0.2288 |
| x2          | +0.1110 | +0.1611 | 23.40 | 0.2388 | 0.1537 |
| x3          | +0.1140 | +0.2125 | 20.77 | 0.1946 | 0.1362 |
| x4          | +0.0306 | +0.0578 | 25.44 | 0.2641 | 0.0935 |
| x3+k+f22    | +0.0646 | +0.1450 | 20.65 | 0.1954 | 0.1126 |

**Conclusion:** total aerosol loading is the first-order control on the twilight DoLP bias —
x4 removes 84% of it (+0.194 -> +0.031) and puts model DoLP at the measured scale (0.25).
The response is strong but noisy/non-monotone at 48 directions (x2 vs x3 medians flip;
shape RMSE non-monotone) — structural field changes + heavy-tailed MC. The joint
(DoLP, intensity, AoP) optimum sits between x2 and x4; the balanced physics-stacked
candidate x3+King+f22r0.7 removes 67% of the bias with the best AoP of any run (20.65 deg)
at +0.02 intensity-RMSE cost. Scan summary archived:
paper/tables/reproducibility/marseille_depol_scan_summary_2026-08-10.csv (+ x4 and
x3_full comparison CSVs).

**Final decomposition of the +0.19 over-polarization**: aerosol loading (dominant,
INPUT-side — the frozen climatological profile understates the real 2022-08-15 aerosol);
molecular King factor -0.023 (model physics, adopt); aerosol nonsphericity small for DoLP
but cleans AoP (adopt as knob); sampling-tail understatement secondary; chain depth /
estimator choice / profile shape aloft — nil. Task #21 CLOSED as a mechanism hunt.

**Follow-on (new work, needs user/externals):** (1) corroborate the hazy-day hypothesis
with AERONET Marseille/OHP AOD for 2022-08-15 (internet, user-side); (2) if corroborated,
build a day-specific aerosol profile (AOD-matched, boundary-layer weighted rather than
uniform xN), rerun the full 683-direction field with King+f22, and re-freeze as the
production tier alongside the chi-fix (task #23's freeze-once plan); (3) the Marseille-vs-
Koomen bias asymmetry (4-6x under the identical chain) stays open as a possible
reduction-side component worth one dedicated audit.

## 2026-08-10 (evening): AERONET corroboration — hazy-day hypothesis REFUTED

Discovery: the login node has outbound HTTPS; AERONET v3 pulled directly via
`print_web_data_v3` (no user lookup needed). Sites: OHP_OBSERVATOIRE L2.0
(inland, 680 m, 75 km NNE; PI Goloub), Toulon L1.5 (coastal, 50 m, 52 km ESE;
PI Piazzola). Carpentras offline (zero points both levels, 08-14..16).
Interpolation AOD(550)=AOD(500)·(1.1)^(−α), α=440–675 Ångström.

Verdict vs frozen column AOD(550)=0.1379 (scan multiples 0.276/0.414/0.552):
- **Toulon (the coastal Marseille analog), Aug 15: afternoon plateau 0.17–0.19;
  last direct-sun point 17:46:39 UTC = 0.1751; next morning first point
  (Aug 16 05:34) = 0.1735 → the 19:14 UTC observation is tightly bracketed at
  AOD(550) ≈ 0.175 ± 0.02 = ×1.27 frozen.** Nowhere near the ×2–×4 the DoLP
  scan requires.
- OHP inland evening declines to 0.072 by 17:49 UTC (below frozen); the ~0.1
  coastal–inland difference is the marine boundary layer Marseille sits in.
- α ≈ 1.15–1.29 (Toulon 08-15): fine-mode; NO dust event (dust ⇒ α<0.5), so
  no external mandate for strong f22 reduction either.
- Aug 14/16 medians ×1.1–1.4 at both sites: regionally coherent MILD
  enhancement, not a haze event.

Implication for #21/#23: the ×2–×4 loading scans were a compensating knob, not
a physical correction. Licensed input-side stack = chi-fix + King(0.0279) +
×1.27 day-matched profile (+f22 0.7 for AoP, weak physical support): removes
only ~0.05–0.06 of the +0.194 signed bias (scan interpolation: ×1.27 ≈ −0.02
to −0.03; King −0.023; f22 −0.008). **Residual ≈ +0.13–0.14 over-polarization
is a genuine open discrepancy** — model MS physics or reduction-side; the
Marseille-vs-Koomen 4–6× bias asymmetry audit is now the prime lead. The
×3+King+f22r0.7 "balanced candidate" is demoted to sensitivity/knob status.
Declared caveat: AERONET = point column vs ~10^3 km elevated twilight path;
but a ×2–×4-equivalent enhancement would have registered at both stations.

Re-freeze (#23) reframed: honest stack + declared residual, NOT an AOD-forced
match. Paper framing B (regression archaeology + honest residual) strengthened.
Artifacts: paper/tables/reproducibility/aeronet_20220815/ (raw pulls, parser,
aod550_summary.csv, README with verdict + PI acknowledgments).

## 2026-08-11: V2 SINGLE RE-FREEZE EXECUTED (task #23 CLOSED)

User-confirmed choices: BL-weighted profile / f22 EXCLUDED (physics-only) /
end-to-end unattended. Mid-execution discoveries that reshaped the freeze:

1. **At-site AERONET found in the frozen dir itself**: `aeronet_direct_sun.csv`
   is Marseille_ATMO (43.306 N, 5.395 E, 65 m, instrument #944) — data ON the
   day. Fresh V3 L2.0 pull is bit-identical (final calibration confirmed).
   Evening plateau 0.1425–0.171 (14:00–17:06 UTC), last point 0.1425; next
   morning 0.2167. Day target revised 0.175 (Toulon) → **0.155** (at-site
   plateau mean = bracket interpolation to 19:14). Toulon 0.175 = regional
   upper bound.
2. **v1 was ALREADY day-informed**: paper_case_provenance.json records
   target_aod_550 = 0.1379 = 0.75×AERONET(0.14393) + 0.25×OpenMeteo(0.12) —
   anchored on the same station's ~17:00 points. v1→v2 (+12%) is an
   evening-extrapolation refinement of the same data, NOT an archaeology.
   First-pass framing corrected in aeronet_20220815/README.md.
3. **ETXTBSY self-inflicted kill (788177)**: rebuilding build_fix3 while a job
   exec'd it. RULE: pinned freeze binary (`build_v2freeze`, md5s in tier
   README) — dev rebuilds stay in build_fix3. Partial outputs wiped; rerun
   clean (RNG case_id-seeded).

**Profile**: BL scale s=1.379584 below 700 m → column 0.155000; non-aerosol
columns bit-identical (`aeronet_20220815/build_day_profile.py`).

**Runs (788181, pinned binary)**: subset 3.7 min, full 683-dir 56 min.
Comparison-CSV metrics (scorer calibrated EXACT vs 08-10 scan table):
- subset: sgn_med +0.1938→+0.1678, AoP 23.37→18.01 (best ever), shape
  0.1731→0.2599, model_dolp_med 0.3877→0.4683 (BL-shadow selective-extinction
  effect: BL aerosol at SZA 96 is unlit — attenuates bright low-DoLP horizon
  light instead of adding unpolarized scatter; uniform-xN scans lowered DoLP)
- full: sgn_med +0.2068→+0.1736, AoP 16.96, shape 0.0827→0.0940, zenith
  buckets +0.350/+0.143/+0.093 → +0.304/+0.106/+0.076 (fingerprint persists)
- day profile contributes only ~−0.003 beyond King → residual ≈ +0.17 is the
  declared open discrepancy (NOT input-side).

**Gate v2 (788184; hardened metric; 22/23 PASS)**:
- `convergence_low_order_metric` knob implemented (MonteCarloDriver.hpp/.cpp
  parser+metadata, ValidationQA evaluateConvergence; default false
  bit-identical — ctest 5/5 + tiny-capped runtime bitcheck vs build_fix2 with
  ZERO non-timing diffs). Peaks on deterministic first+second order = exactly
  0; flux stays total-field.
- **Rozenberg 0.0851442 ≈ claim 0.0849** (true h-6 + chifix + King in-gate).
- **KOOMEN FULL PASS, FIRST TIME**: RMSE 0.040091 (was ~0.10),
  median_dolp_abs 0.0160022 < 0.02 (was 0.0515 chifix / 0.0309 legacy). Same
  h-6 input archaeology applied coherently to the same scenario family —
  Koomen claim-era numbers evidently also from the unrecorded h-6 profile.
  The 4-6x Marseille-vs-Koomen asymmetry lead sharpens: identical physics
  stack passes Koomen and leaves Marseille +0.17.
- **Declared deviation (the 1 fail)**: convergence_flux_rel 0.0247044 vs 0.02,
  BUDGET-INDEPENDENT (0.0239651 at 256 photons /788183/, 0.0247044 at 1024
  /788184/) — heavy-tailed higher-order estimator; legacy 0.006 was estimator
  suppression, not better convergence. Photon budget kept at 1024.

**Frozen artifacts**: tier `data/paper_cases/frozen_..._v2/` (profile, README
with binary md5s, configs/); archives in `paper/tables/reproducibility/`:
marseille_v2_freeze_summary_2026-08-11.csv, gate report json, both comparison
CSVs, aeronet_20220815/. Docs updated: STATUS_README (2026-08-11 addendum +
supersession note), monte_carlo_cpp/README, MODEL_ASSUMPTIONS. New tool:
tools/score_measurement_comparison.py.

**Open after freeze**: (a) Marseille reduction-chain audit (prime lead for the
+0.17 residual); (b) paper framing A/B decision (B strengthened: honest
residual + Koomen vindication); (c) optional f22/aerox sensitivity table in
paper from 08-10 scans.
