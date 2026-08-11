# GW Si8 Optimization Step Log

## Scope
- Work directory: `GW/test_profile/test_nok_si8`
- Goal: explain and improve server-side GW runtime, especially `input` phase bottlenecks.
- Rule: every optimization step must be logged here with reproducible details.

## Baseline (before changes)
- Date: 2026-05-06
- Source: MATLAB profiler outputs `profile_input/file0.html` and `profile_qp/file0.html`.
- Observed wall-time split (server):
  - `input_driver + load(config)`: ~162 s (dominant)
  - `qp_cohsex`: ~6 s
- Key hotspots inside `input`:
  - `adaptiveisdf`: ~121 s total
  - `adaptive_weight`: ~64 s total (~44 s self)
  - `isdf_schur_update`: ~52 s total (~49 s self)
  - `prod`: ~22.6 s total
- Preliminary conclusion:
  - Fast BLAS kernel alone cannot strongly reduce end-to-end runtime.
  - Main bottleneck is adaptive ISDF workflow and many high-frequency small calls.

## Step Log

### Step 0 - Establish profiling baseline
- Status: done
- What changed:
  - No code changes yet.
  - Collected and interpreted hotspot distribution.
- Why:
  - Need quantitative bottleneck map before optimization.
- Evidence:
  - `profile_input/file0.html`
  - `profile_qp/file0.html`
- Next:
  - Step 1 will add large BLAS baseline for `gemm/trsm/trmm`.

### Step 1 - Add large BLAS benchmark (`gemm/trsm/trmm`)
- Status: done
- What changed:
  - Added a new MATLAB benchmark script for large matrices:
    - `gemm`: `A * B` with square `n x n`
    - `trmm`: `U * X` with upper-triangular `U`
    - `trsm`: `U \ X` with upper-triangular `U`
  - Script includes warmup, repeated timings, `best/mean` stats, and GFLOPS summary table.
  - Default sweep:
    - dtype: `single`, `double`
    - `n`: `[1024, 2048, 4096]`
    - `rhs`: `[256, 512, 1024]` (for triangular ops)
- Why:
  - Build a dedicated large-matrix BLAS baseline before optimizing GW hotspots.
- Files touched:
  - `GW/test_profile/test_blas/test_large_blas_gflops.m`
- Validation:
  - Static code review completed.
  - Runtime benchmark not executed yet in this step.
- Result:
  - Benchmark harness created and ready for server/laptop side-by-side runs.
- Notes/Risks:
  - `n=4096` in `double` may consume substantial memory; reduce `n_list` if needed.
- Next:
  - Run script on server first, collect summary table, then run on laptop with same settings.

### Step 2 - Summarize small-vs-large GFLOPS loss
- Status: done
- What changed:
  - Consolidated benchmark interpretation into a dedicated markdown note under `test_blas`.
  - Quantified GFLOPS loss of small/narrow workload relative to large GEMM peak.
- Why:
  - Need a clear, reusable answer to "small scale loses how much GFLOPS vs large scale".
- Files touched:
  - `GW/test_profile/test_blas/blas_efficiency_summary.md`
- Validation:
  - Used server terminal benchmark outputs from:
    - `test_isdf_prod_batch_gflops`
    - `test_large_blas_gflops`
  - Computed loss as `(1 - G_small / G_ref)`.
- Result:
  - Small/narrow throughput generally loses about `~82% to ~98%` vs large GEMM peak.
- Notes/Risks:
  - Current summary is server-side only; laptop-side same script run is still needed for strict apples-to-apples comparison.
- Next:
  - Run `test_large_blas_gflops` and `test_isdf_prod_batch_gflops` on laptop with identical settings and append cross-machine comparison.

### Step 3 - Add `Nold`-wise batch recommendation table
- Status: done
- What changed:
  - Extended BLAS summary with an explicit `Nold -> best batch` table.
  - Added a practical piecewise starting policy and autotune suggestion for `adaptive_weight`.
- Why:
  - Need an actionable batch choice rule, not only throughput observations.
- Files touched:
  - `GW/test_profile/test_blas/blas_efficiency_summary.md`
- Validation:
  - Derived directly from current server benchmark table (`test_isdf_prod_batch_gflops`).
  - Cross-checked best-batch choices for all four `Nold` points.
- Result:
  - Current best-by-`Nold` mapping:
    - `1e4 -> 256`, `5e4 -> 128`, `1e5 -> 64`, `2e5 -> 256`.
- Notes/Risks:
  - Mapping may shift with CPU model, thread settings, and dtype.
- Next:
  - Implement batch-selection logic (range-based or cached autotune) inside `adaptive_weight`.

### Step 4 - Fix `data >2GB` save warning in `input_driver`
- Status: done
- What changed:
  - Removed unintended bare `save()` call in `driver_profile/input_driver.m`.
- Why:
  - `save()` with no arguments writes `matlab.mat` and attempts to dump full workspace.
  - Large `data` variable triggered warning:
    - `Variable 'data' was not saved. For variables larger than 2GB use MAT-file version 7.3 or later.`
- Files touched:
  - `GW/driver_profile/input_driver.m`
- Validation:
  - Static verification: warning stack in log points to `input_driver (line 28)`, exactly matching removed call.
  - Existing explicit saves already handle large data safely where needed:
    - `save(fNamedata, 'data', '-v7.3', '-nocompression')`
- Result:
  - The spurious `matlab.mat` workspace save path is removed; this warning should disappear in next run.
- Notes/Risks:
  - None for normal workflow; this only removes debug-like side effect.
- Next:
  - Re-run `test_nok_si32` and confirm warning is gone from `log/test_nok_si32.log`.

## Per-step Template
Use this block for each new step:

```md
### Step N - <short title>
- Status: planned | in_progress | done
- What changed:
  - <code/config changes>
- Why:
  - <motivation>
- Files touched:
  - `<path>`
- Validation:
  - <how validated, timing/profiler diffs>
- Result:
  - <before vs after numbers>
- Notes/Risks:
  - <known tradeoffs>
- Next:
  - <next action>
```
