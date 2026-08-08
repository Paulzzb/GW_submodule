# Si8_sc — small supercell SC_ISDF demo

2×2×1 supercell over `Si2_uc`, abstracted from `test_profile/test_SC/Si8_from_Si2` + `test_formal/Si8_formal`.

| Item | Choice |
|------|--------|
| Groundstate | `formal` from `../Si2_uc/qe.save` xml (no large SC QE hdf5) |
| ISDF | `use_sc_isdf` + `sc_adaptive`, source `../Si2_uc/SAVE` |
| Frequency | `-2` (Gamma COHSEX path required by formal) |
| Size | FFT ~36×36×18, 32 bands — keep this demo light |

## Prerequisite

Run `Si2_uc` first so `../Si2_uc/SAVE` exists (or use `example/run_sc`).

## Run

```matlab
cd example
run_sc('sc')     % SC only (UC SAVE must exist)
run_sc           % UC then SC + qp.launcher
```
