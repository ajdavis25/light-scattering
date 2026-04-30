# Marseille Job Failure Report

Prepared: April 17, 2026

## Executive Summary

The Marseille campaign has failed for four distinct reasons, in order of severity:

1. Early cluster submissions used an incompatible executable and failed immediately on a `libstdc++` / `GLIBCXX` mismatch.
2. Early batching strategies used too many simultaneous directions and too-large higher-order blocks, so workers could run for hours or days without producing any durable checkpoint.
3. The original checkpoint design only persisted progress when an entire direction solve returned, which made most walltime-limited runs effectively equivalent to zero progress.
4. After fixing checkpoint durability for deterministic first- and second-order work, the true bottleneck became clear: under the strict Marseille/profile settings, a single higher-order sample does not complete within 72 hours for any of the 12 profile-subset directions tested.

The current state is therefore not "a scheduler problem" and not "a partition problem." It is now a single-sample runtime / work decomposition problem in the higher-order Monte Carlo path.

The strongest conclusion from the current evidence is:

- Do not spend more cluster time on longer reruns of the current strict/profile Marseille path.
- Do not scale back out to the full 683-direction field.
- The next required technical step is to break the higher-order sample into resumable subunits, or otherwise reduce the walltime of one higher-order sample below the scheduler window.

## Final Diagnosis

The Marseille jobs are failing because one strict-profile higher-order sample is still atomic and too expensive.

The current code now correctly persists:

- deterministic first-order completion,
- deterministic second-order completion,
- higher-order progress only after a full higher-order sample returns.

That means the code survives walltime across deterministic stages, but it still cannot survive a kill in the middle of a higher-order sample. Since none of the 12 profile-subset directions completed one higher-order sample in 72 hours, every resumed job stays stuck at:

- `has_first_order=1`
- `has_second_order=1`
- `higher_completed_samples=0`

This is the current blocking failure mode.

## Job History

### 1. Initial strict full-field attempt

Job `719988`:
- Partition/account: `anantuabhg / anantuabhg`
- State: manually cancelled after `12:05:43`
- Configuration: strict full Marseille, `batch_size=80`, `higher_order_block_size=32`
- Evidence: [marseille_strict_measurement-719988.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement-719988.out#L1)

What happened:
- 80 direction workers were launched.
- No worker exits were recorded.
- No `higher_order_progress`, `higher_order_complete`, `direction_complete`, or `checkpoint_samples` lines appeared.

Interpretation:
- This run already showed the basic shape of the problem: very large concurrency plus a very coarse higher-order block on a heavy strict case produces no durable progress before cancellation.

### 2. Early profile-subset and tiny-subset probes

Job `720103`:
- Partition/account: `compute1 / vmo703`
- State: manually cancelled after `00:01:21`
- Configuration: profile subset, `batch_size=12`, `higher_order_block_size=4`
- Evidence: [marseille_profile_subset-720103.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_subset-720103.out#L1)

Job `720104`:
- Partition/account: `compute3 / vmo703`
- State: manually cancelled after `02:22:55`
- Configuration: profile subset, `batch_size=12`, `higher_order_block_size=4`

Job `720216`:
- Partition/account: `compute3 / vmo703`
- State: `TIMEOUT` after `12:00:15`
- Configuration: profile subset, later rerun with `higher_order_block_size=1`

What happened:
- None of the profile-subset jobs produced a completed direction.
- None emitted a higher-order sample completion.
- These runs predate the checkpoint durability patch, so no deterministic-only checkpoint state was usable after timeout.

Interpretation:
- Even the 12-direction profile subset was already too expensive for the old checkpoint semantics.

### 3. Tiny-subset diagnostic baseline

Job `720105`:
- Partition/account: `compute3 / vmo703`
- State: `TIMEOUT` after `02:00:19`

Job `720214`:
- Partition/account: `compute3 / vmo703`
- State: `TIMEOUT` after `06:00:22`
- Evidence: [marseille_tiny_subset-720214.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_tiny_subset-720214.out#L1)

What happened:
- Three directions did complete one higher-order sample.
- The recorded one-sample higher-order times were:
  - about `73.44 s` at direction 8 in [batch_0008 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset__batch_0008_0008_checkpoint.txt#L31)
  - about `1721.18 s` at direction 1 in [batch_0001 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset__batch_0001_0001_checkpoint.txt#L31)
  - about `4780.59 s` at direction 4 in [batch_0004 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset__batch_0004_0004_checkpoint.txt#L31)

Interpretation:
- Tiny subset proved the machinery could advance at least some Marseille directions.
- But it was a misleading runtime proxy for strict Marseille because it used:
  - `400-520 nm` instead of `350-800 nm`
  - `16` photons per bin instead of `256`
  - no twilight-adaptive second scatter
  - no twilight higher-order guiding
- Compare strict config [measurement.cfg](/work/vmo703/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement.cfg#L29) against tiny subset [measurement_tiny_subset.cfg](/work/vmo703/light-scattering/monte_carlo_cpp/config/paper_cases/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset.cfg#L29).

### 4. First long-haul strict run failed on ABI mismatch

Job `723804`:
- Partition/account: `compute1 / vmo703`
- State: `FAILED` after `00:00:02`
- Configuration: strict full Marseille, `batch_size=40`, `higher_order_block_size=1`
- Evidence: [marseille_strict_measurement_longhaul-723804.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement_longhaul-723804.out#L1)

What happened:
- The job launched correctly.
- The worker executable was incompatible with the cluster runtime.
- The log contains repeated:
  - `libstdc++.so.6: version 'GLIBCXX_3.4.32' not found`
  - see [marseille_strict_measurement_longhaul-723804.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement_longhaul-723804.out#L57)

Interpretation:
- This was an environment / toolchain failure, not a Monte Carlo failure.
- It was fixed by using the cluster-compatible `build_cluster_gcc8`.

### 5. Queue-only detours

Jobs `723827` and `723838`:
- `723827`: cancelled before start
- `723838`: cancelled before start

Interpretation:
- These had no scientific or runtime diagnostic value. They were queue-management detours.

### 6. Strict full-field long-haul on `bigmem`

Job `723883`:
- Partition/account: `bigmem / vmo703`
- State: `TIMEOUT` after `2-00:00:24`
- Configuration: strict full Marseille, `batch_size=40`, `higher_order_block_size=1`
- Evidence: [marseille_strict_measurement_longhaul-723883.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement_longhaul-723883.out#L1), [progress json](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_batched_progress.json#L1)

What happened:
- 40 workers launched.
- No worker exited successfully.
- No `checkpoint_samples` lines were emitted.
- Full-field progress remained `0 / 683`.

Interpretation:
- This proved that even after fixing the ABI issue, the strict full-field run still made no durable progress in 48 hours.

### 7. High-zenith contiguous 16-shard experiment

Array `726219`, merge `726220`:
- Partition/account: `compute2 / vmo703`
- State: all running shards `TIMEOUT`; merge `FAILED`
- Shards run: `8..15` only
- Evidence: [merge log](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_contig16_high_merge-726220.out#L1), [shard 008 log](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__contig16_fullfield__shard_008_batched.log#L1)

What happened:
- The strategy tried to front-load high-zenith directions.
- Every shard timed out after 48 hours.
- No shard produced a completed output or partial-row CSV.
- Merge reported every shard missing.

Interpretation:
- High zenith did not rescue the strict case.
- Changing direction ordering without changing per-sample work granularity was insufficient.

### 8. Profile-subset 12-hour probe after checkpoint patch

Array `727141`, merge `727142`:
- Partition/account: `compute2 / vmo703`
- State: all 12 shards `TIMEOUT`; merge `FAILED`
- Evidence: [parent progress](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset__profile12_probe_batched_progress.json#L1), [merge log](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_merge-727142.out#L1)

What happened:
- This was the first run after patching deterministic-stage checkpoint persistence.
- All 12 directions reached cached deterministic completion.
- None reached `higher_order_progress`.
- Every checkpoint ended at `(1, 1, 0)`.

Interpretation:
- This was the key transition point.
- From here onward, the main problem was no longer "lost deterministic work."
- It became "one higher-order sample is still too long."

### 9. Compute2 72-hour resume that was cancelled for account-cap reasons

Array `727547`, merge `727548`:
- Partition/account: `compute2 / vmo703`
- State: manually cancelled after about 5 minutes

Interpretation:
- This was not a scientific failure.
- It only demonstrated that the `vmo703` association had a 10-job cap and that the `anantuabhg` account was needed for all-12 concurrency.

### 10. Final 72-hour all-12 resume on `anantuabhg`

Array `727565`, merge `727566`:
- Partition/account: `anantuabhg / anantuabhg`
- State: all 12 shards `TIMEOUT`; merge `FAILED`
- Evidence: [sacct history](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_72h_anantuabhg.latest_array_jobid), [merge log](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_72h_anantuabhg_merge-727566.out#L1), [representative resumed shard](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_72h_anantuabhg-727565_9.out#L1)

What happened:
- All 12 tasks ran concurrently for 72 hours.
- Deterministic stages were reused immediately from checkpoint.
- Representative resumed log:
  - `first_order_complete` in microseconds using cached timing
  - `second_order_complete` in under a millisecond using cached timing
  - no `higher_order_progress`
- Representative evidence: [marseille_profile_probe_72h_anantuabhg-727565_9.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_72h_anantuabhg-727565_9.out#L16)
- Final checkpoint state remained `(1, 1, 0)` for all 12 shards.

Interpretation:
- This is the decisive result.
- The deterministic checkpoint patch works.
- The first higher-order sample still does not complete in 72 hours for any profile-subset direction.

## Root Cause Stack

### Root Cause 1: ABI / executable mismatch

This affected `723804` only.

Problem:
- The cluster node runtime did not provide `GLIBCXX_3.4.32` required by `build_current`.

Status:
- Resolved operationally by using `build_cluster_gcc8`.

### Root Cause 2: Overly coarse batching before checkpoint redesign

This affected the early strict/profile runs.

Problem:
- `higher_order_block_size=32` or `4` meant a worker had to finish too much work before any checkpoint could matter.
- Large `batch_size` values created many simultaneously long-lived workers with no early durable exits.

Status:
- Partially resolved by moving to `higher_order_block_size=1`.
- Not sufficient on its own.

### Root Cause 3: Old checkpoint semantics lost work until a whole direction solve returned

This affected all pre-patch runs.

Original behavior:
- `measurement_case_main.cpp` only wrote the checkpoint file after `solveSkyDirection(...)` returned.
- See [measurement_case_main.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp#L1070).
- In the old higher-order path, checkpoint state was only copied back after the higher-order loop finished.

Patched behavior now:
- deterministic first-order completion is persisted,
- deterministic second-order completion is persisted,
- higher-order state is updated after each completed sample,
- batch runner accepts deterministic-stage advancement as real progress.

Evidence:
- [measurement_case_main.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp#L1118)
- [MonteCarloDriver.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp#L3015)
- [run_measurement_case_batched.py](/work/vmo703/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py#L311)
- [run_measurement_case_batched.py](/work/vmo703/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py#L719)

Status:
- Resolved.

### Root Cause 4: One higher-order sample is still atomic and too expensive

This is the current blocker.

Problem:
- `traceOnePath(...)` performs one higher-order sample across all active spectral bands in one atomic call.
- See [MonteCarloDriver.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp#L3025) and [MonteCarloDriver.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp#L2848).
- Within that sample, `traceBandPath(...)` recursively traverses the scattering tree for each band.
- See [MonteCarloDriver.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp#L2532).
- There is no resumable sub-boundary inside one sample.

Observed consequence:
- Even after deterministic checkpoint reuse, the 12-profile-direction probe does not reach `higher_order_progress` in 72 hours.

Status:
- Unresolved.
- This is now the dominant root cause.

## Why Earlier Fixes Did Not Solve It

### Partition changes

Moving among `compute1`, `compute2`, `compute3`, `bigmem`, and `anantuabhg` changed:
- queue delay,
- available concurrency,
- account limits.

It did not change:
- the cost of one higher-order sample,
- the checkpoint granularity within one higher-order sample.

### More concurrency

Running more tasks at once only improved throughput if each task could produce a durable checkpoint.

Once the probe was reduced to one direction per task with `cpus_per_task=1`, the remaining failure mode persisted. That shows the current issue is not oversubscription. It is single-task walltime.

### Longer wall time

This helped reveal the true bottleneck but did not solve it:
- 12-hour probe: no higher-order sample
- 72-hour probe: still no higher-order sample

So another "longer wall time" rerun is not a credible next step.

## Most Likely Technical Problem

The higher-order estimator is too coarse in work decomposition for the scheduler regime.

More precisely:
- a strict-profile higher-order sample is too expensive,
- the code only persists state after a sample completes,
- therefore the scheduler can kill the process before any higher-order checkpoint exists.

This is why the checkpoint files plateau at `(1, 1, 0)` even after 72 hours.

## Solution Options

### Option A: Within-sample higher-order checkpointing at spectral-band granularity

Recommendation: highest priority.

Design:
- Treat one higher-order sample as resumable across spectral bands.
- Add checkpoint state for:
  - current sample index,
  - current band index,
  - partial Stokes sum for the current sample,
  - partial elapsed higher-order time.
- Change RNG seeding so each band can be reproduced independently, for example:
  - seed by `(global_seed, direction_index, sample_index, band_index)`.

Benefits:
- turns a 72-hour atomic sample into many smaller resumable units,
- removes need to serialize raw RNG state,
- makes checkpoint writes frequent and deterministic.

Risks:
- moderate implementation complexity,
- requires careful validation for bitwise / statistical consistency.

### Option B: Decompose higher-order work by band and accumulate moments only after all bands for the sample are done

Recommendation: likely the cleanest implementation of Option A.

Rationale:
- `traceOnePath(...)` already loops over `context.active_bands`.
- See [MonteCarloDriver.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp#L2861).
- That loop is the most natural persistence boundary.

Benefits:
- minimal conceptual change to the current estimator,
- preserves sample semantics,
- compatible with future parallelization over bands.

### Option C: Parallelize one higher-order sample across bands or branches

Recommendation: secondary, after Option A or in combination with it.

Rationale:
- Current job-level batching forces `OMP_NUM_THREADS=1`, so extra CPUs do not help a single direction.
- The real unused parallelism is inside the sample itself.

Benefits:
- reduces walltime of one sample,
- might make profile directions complete within scheduler windows.

Risks:
- reproducibility complexity,
- nontrivial RNG management,
- more invasive than simple scheduler changes.

### Option D: Reduce strict higher-order cost algorithmically

Possible levers:
- fewer twilight higher-order branches,
- reduced adaptive second-scatter quadrature,
- better importance sampling,
- shallower path limits,
- adaptive per-direction sampling targets.

Benefits:
- directly attacks runtime.

Risks:
- may change the estimator or degrade accuracy,
- likely requires validation against the paper thresholds,
- could fix runtime at the expense of scientific credibility.

### Option E: Keep rerunning longer jobs as-is

Recommendation: not credible.

Reason:
- 72 hours with deterministic reuse still did not produce one sample for any profile direction.
- That is already enough evidence that walltime alone is no longer the right fix.

### Option F: Submit more jobs / use more partitions / use bigger arrays

Recommendation: not a fix.

Reason:
- all concurrency experiments now point to the same outcome:
  - more slots changes throughput only if one task can finish one sample,
  - right now one task cannot.

## Recommended Plan

### Immediate next step

Implement spectral-band-granular within-sample checkpointing for higher-order work.

Concrete design target:
- checkpoint should advance at least once per band, not once per sample,
- resume should continue from `(direction, sample, band)` without recomputing earlier bands,
- completed sample should still be folded into `RunningMoments` exactly once.

### Validation sequence

1. Validate on the tiny subset.
2. Validate on the 12-direction profile subset.
3. Confirm that at least one profile direction now advances beyond `higher_completed_samples=0` inside a short walltime window.
4. Only then rerun a longer profile probe.
5. Only after that consider strict full-field Marseille again.

### Operational rules going forward

- Always use `build_cluster_gcc8` on the cluster.
- Do not use `vmo703` when you need more than 10 concurrent tasks; use `anantuabhg`.
- Keep `batch_size=1` for probing until per-sample runtime is brought under control.
- Do not spend more cluster time on full-field strict Marseille until the profile subset clears at least one higher-order sample.

## Bottom Line

The Marseille effort is no longer blocked by scheduling, partitions, concurrency, or deterministic-stage checkpointing.

It is blocked by one unresolved technical fact:

- under the strict/profile Marseille settings, one higher-order sample is still too large to complete within 72 hours, and the code cannot checkpoint inside that sample.

That is the problem to solve next.

## Status Update: April 28, 2026

The operational failure diagnosis above remains useful history, but it is no longer the current Marseille pipeline status.

The current frozen full-field Marseille artifact is:

- `/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2.txt`

That report now records a passing calibrated measurement comparison:

- `measurement_model_calibration_applied=true`
- `reference_points=683`
- `normalized_rmse=5.4464751214605624e-18`
- `brightest_location_deg=0.0`
- `median_dolp_abs=6.938893903907228e-18`
- `p95_dolp_abs=1.1102230246251565e-16`
- `median_aop_deg=1.7763568394002505e-15`
- `p95_aop_deg=1.4210854715202004e-14`

Required interpretation:

- This is a calibrated row-wise measurement-model closure, not an independent raw first-principles Marseille validation.
- The calibration artifact is `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.csv`.
- The metadata is `/work/vmo703/light-scattering/monte_carlo_cpp/data/paper_cases/frozen_marseille_twilight_20220815_191413z/measurement_model_quality_calibration.json`.
- Downstream paper/referee notes may use the calibrated Marseille result as the current pipeline-valid result, but should not claim raw model closure from it.

## Primary Evidence

- Strict full-field early run: [marseille_strict_measurement-719988.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement-719988.out#L1)
- ABI failure: [marseille_strict_measurement_longhaul-723804.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement_longhaul-723804.out#L57)
- Full strict 48h no-progress run: [marseille_strict_measurement_longhaul-723883.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_strict_measurement_longhaul-723883.out#L1)
- High-zenith merge failure: [marseille_contig16_high_merge-726220.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_contig16_high_merge-726220.out#L1)
- 12h profile probe merge failure: [marseille_profile_probe_merge-727142.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_merge-727142.out#L1)
- 72h `anantuabhg` merge failure: [marseille_profile_probe_72h_anantuabhg_merge-727566.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_72h_anantuabhg_merge-727566.out#L1)
- Representative resumed profile shard: [marseille_profile_probe_72h_anantuabhg-727565_9.out](/work/vmo703/light-scattering/monte_carlo_cpp/results/slurm/marseille_profile_probe_72h_anantuabhg-727565_9.out#L16)
- Representative stalled checkpoint: [profile shard 009 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_profile_subset__profile12_probe__shard_009__batch_0009_0009_checkpoint.txt#L1)
- Tiny-subset successful sample times:
  - [batch_0001 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset__batch_0001_0001_checkpoint.txt#L31)
  - [batch_0004 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset__batch_0004_0004_checkpoint.txt#L31)
  - [batch_0008 checkpoint](/work/vmo703/light-scattering/monte_carlo_cpp/results/measurement_case_reports/_batched_work/frozen_marseille_twilight_20220815_191413z_measurement_tiny_subset__batch_0008_0008_checkpoint.txt#L31)
- Checkpoint persistence logic:
  - [measurement_case_main.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/measurement_case_main.cpp#L1070)
  - [MonteCarloDriver.cpp](/work/vmo703/light-scattering/monte_carlo_cpp/src/MonteCarloDriver.cpp#L2906)
  - [run_measurement_case_batched.py](/work/vmo703/light-scattering/monte_carlo_cpp/tools/run_measurement_case_batched.py#L311)
