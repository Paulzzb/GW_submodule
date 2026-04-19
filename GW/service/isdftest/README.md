# ISDF Test Sandbox

This folder is an isolated sandbox for testing ISDF index-selection logic.

## Purpose
- Validate selection rules without changing production flow.
- Focus first on coarse-grid subset selection from FFT R-grid.

## Suggested workflow
1. Restore service stage data.
2. Read `fft_data`, `isdf_data`, and `pair_symmetry` objects from managers.
3. Run test entry script and inspect printed diagnostics.
4. Compare selected index count against expected threshold logic.

## Files
- `run_coarse_selection_test.m`: test entry for coarse index-selection checks.
