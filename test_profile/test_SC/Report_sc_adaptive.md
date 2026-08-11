# test_SC Validation Report — sc_adaptive (Si + LiH)

Generated: 2026-06-11 19:51:50

Plain `SC_ISDF` results are in [Report.md](./Report.md).

## Configuration

- `compute_vc = .true.` and `compute_vn = .true.` in Si2 / Si8 / Si16 `test` files
- `energy_band_index_max` = 2 × occupied bands (Si2: 8, Si8: 32, Si16: 64)
- **Supercell cases use `sc_adaptive = .true.`** in `&SUPERCELL`: SC-replicated unit-cell grid as adaptive seed, then `adaptiveisdf` refinement (`adaptive_double` backend)
- LiH supercell cases: `compute_vn = .true.` only; LiH_666 stopped / not included

## HF Validation Summary (sc_adaptive)

| Case | Type | id | nisdf | scheme | ‖ΔE‖₂/‖E_HF‖₂ |
|------|------|-----|-------|--------|----------------|
| Si2 | vc | 2 | 19 | adaptive | 6.42760221e-07 |
| Si2 | vn | 4 | 26 | adaptive | 2.16195033e-06 |
| Si8_from_Si2 | vc | 2 | 76 | adaptive | 5.43428518e-01 |
| Si8_from_Si2 | vn | 4 | 104 | adaptive | 1.00505878e+00 |
| Si16_from_Si2 | vc | 2 | 152 | adaptive | 5.67230507e-01 |
| Si16_from_Si2 | vn | 4 | 208 | adaptive | 3.97002021e-01 |
| LiH_222_from_111 | vn | 2 | 808 | adaptive | 1.88715614e-01 |
| LiH_444_from_111 | vn | 2 | 6464 | adaptive | 4.07796074e-01 |

## Observations

- **Si2 unit cell:** reference adaptive runs (no supercell); machine-precision HF agreement (relative L2 ≲ 10⁻⁶).
- **Si8 / Si16 SC + sc_adaptive:** HF errors **match** plain `SC_ISDF` (`adaptive_sc`) runs on relative L2; adaptive refinement adds **0** centroids because the replicated grid already saturates the ISDF-ratio cap / meets the loss threshold.
- **LiH_222 / LiH_444:** same behaviour — full replicated seed, 0 adaptive adds, HF accuracy unchanged vs plain SC.
- **Implementation:** `SC_ISDF_prepare_adaptive_seed` deduplicates `fine_grid_lin`, sets `interp_scheme='coarse'`, and hands off to `adaptiveisdf`. A small Cholesky jitter in `isdf_schur_update('init')` handles near-singular CCH from redundant SC replicas without altering the Schur update logic.

## sc_adaptive driver flow

```
SC_ISDF → set_nrange → SC_ISDF_prepare_adaptive_seed → adaptive_double.adaptiveisdf → gen_tildeVq → validate_hf
```

Input flag: `&SUPERCELL sc_adaptive = .true.` (requires `use_sc_isdf = .true.`).

## Validation reports on disk

| Case | Report file |
|------|-------------|
| Si8_from_Si2 vc | `Si8_from_Si2/isdf_validate_HF_id2.txt` |
| Si8_from_Si2 vn | `Si8_from_Si2/isdf_validate_HF_id4.txt` |
| Si16_from_Si2 vc | `Si16_from_Si2/isdf_validate_HF_id2.txt` |
| Si16_from_Si2 vn | `Si16_from_Si2/isdf_validate_HF_id4.txt` |
| LiH_222_from_111 vn | `LiH_222_from_111/isdf_validate_HF_id2.txt` |
| LiH_444_from_111 vn | `LiH_444_from_111/isdf_validate_HF_id2.txt` |
