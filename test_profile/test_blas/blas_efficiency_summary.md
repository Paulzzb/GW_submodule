# BLAS Efficiency Summary (test_blas)

Date: 2026-05-06  
Scope: compare "small/narrow workload" vs "large workload" in GFLOPS on server.

## Data Sources
- Small/narrow proxy:
  - `test_isdf_prod_batch_gflops` summary (from terminal output).
  - Kernel: `(conj(Psi)*conj(psi)') .* (Phi*phi')`, with `k1=k2=64`, `batch in [16,64,128,256]`, `Nold in [1e4,5e4,1e5,2e5]`.
- Large workload baseline:
  - `test_large_blas_gflops` summary (from terminal output).
  - Main reference peak:
    - `single GEMM n=4096`: **2199.8 GFLOPS (best)**.
  - Extra references:
    - `single TRMM n=4096,rhs=1024`: **884.6 GFLOPS (best)**.
    - `single TRSM n=4096,rhs=1024`: **611.8 GFLOPS (best)**.

## What "loss" means
For a small-case throughput `G_small` and large-case reference `G_ref`:

`loss(%) = (1 - G_small / G_ref) * 100%`

This section uses `G_ref = 2199.8` (single large GEMM peak).

## Small vs Large: GFLOPS Loss (key points)

### Best small-case points (from `test_isdf_prod_batch_gflops`)
- Best observed small-case throughput: **401.1 GFLOPS** (`Nold=50000,batch=128`)
  - vs large GEMM peak 2199.8: **81.8% loss** (only 18.2% of peak).
- Next strong small point: **350.6 GFLOPS** (`Nold=100000,batch=64`)
  - vs large GEMM peak: **84.1% loss**.
- Many points around **90~132 GFLOPS**
  - vs large GEMM peak: roughly **94.0% ~ 95.9% loss**.
- Worst points around **36~48 GFLOPS**
  - vs large GEMM peak: roughly **97.8% ~ 98.4% loss**.

### By Nold (best over batch)
- `Nold=10000`: best 340.5 GFLOPS -> **84.5% loss**
- `Nold=50000`: best 401.1 GFLOPS -> **81.8% loss** (best among all Nold)
- `Nold=100000`: best 350.6 GFLOPS -> **84.1% loss**
- `Nold=200000`: best 131.1 GFLOPS -> **94.0% loss**

## Batch selection guide (from current server data)

### Best `batch` by `Nold` (using `gflops_best`)
| Nold | best batch | best GFLOPS | note |
|---:|---:|---:|---|
| 10,000 | 256 | 340.5 | Larger batch helps at this size |
| 50,000 | 128 | 401.1 | Global best in this sweep |
| 100,000 | 64 | 350.6 | Mid-size batch wins |
| 200,000 | 256 | 131.1 | Large Nold penalizes all batches; 256 still best |

### Recommended starting policy for `adaptive_weight`
- Start with a piecewise default:
  - `Nold <= 2e4`: try `batch=256`
  - `2e4 < Nold <= 7e4`: try `batch=128`
  - `7e4 < Nold <= 1.5e5`: try `batch=64`
  - `Nold > 1.5e5`: try `batch=256` (from current data, but verify)
- Add a lightweight warmup autotune on first call:
  - test candidate batches `[64, 128, 256]` on a small sample,
  - cache winner by `(Nold range, dtype, k1, k2)`,
  - reuse cached batch to avoid repeated tuning overhead.

### Important caveat
- The optimal batch is **not monotonic** with `Nold` in this dataset.
- Do not hard-code a single global batch value; use range-based or autotuned selection.

### Interpretation
- Relative to large dense GEMM, narrow/small-style workload is typically losing **~82% to ~98%** GFLOPS.
- Even the best small case is only about **1/5.5** of large GEMM peak.
- For many practical points, effective throughput is closer to **1/15 to 1/45** of large GEMM peak.

## Why this is expected
- Small/narrow kernels have lower arithmetic intensity and poorer cache reuse.
- High-frequency small calls amplify MATLAB/runtime dispatch and memory traffic overhead.
- Mixed operations (two narrow GEMMs + elementwise product + data movement) prevent sustained peak BLAS throughput.
- Therefore, machine-level peak (seen in large GEMM) does not translate to end-to-end GW input stage.

## Practical takeaway for GW optimization
- Do not use large GEMM peak alone to estimate total speedup.
- Priority should be:
  1. reduce small-call count,
  2. increase batching/block size where possible,
  3. merge kernels to reduce intermediate memory traffic,
  4. keep thread/NUMA settings stable during profiling.
