# test_Si_isdftest — +isdf vs +isdftest driver comparison

Copy of `test_Si` (without `profile_results` and old adaptive sweep outputs) for comparing
`+isdf` and `+isdftest` on the same Si QE dataset.

## Quick run

From MATLAB, in this directory:

```matlab
run_isdf_vs_isdftest
```

Force rebuild of `SAVE/` from `Si.save/` (skip cached GWinput):

```matlab
run_isdf_vs_isdftest(true)
```

## What it does

1. **Phase 1 — `input_driver('./test')`**
   - Loads/builds `SAVE/`, runs `service_driver` → **`isdf.driver`**
   - Moves `isdf_validate_HF_id*.txt` and `adaptiveisdf_id*.txt` to **`results/isdf/`**

2. **Phase 2 — `isdftest.driver`**
   - Restores relay stage from `test_relay_stage.mat`
   - Runs **`isdftest.driver`** with `pwd = results/isdftest/` so HF reports land there

3. **`results/comparison_summary.txt`**
   - Lists both report sets and extracts `sum(E_HF)` / `sum(E_HF_ISDF)` for side-by-side check

## Prerequisites

- `test`, `Si.save/`, MATLAB path via `QPstartup`
- `isisdf = 1` in `test` (already set)

## Notes

- `+isdf` and `+isdftest` use **separate** manager pools; phase 2 does not overwrite phase 1 data.
- Phase 1 may take several minutes (adaptive ISDF + HF validation for vc/vn/nn).
