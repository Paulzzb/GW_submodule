# Si2_uc — unit-cell ISDF source for SC demo

Small Si primitive cell (2 atoms, FFT 18³, 3 k-points).  
Abstracted from `test_profile/test_SC/Si2`, using the same groundstate family as `Si_k/qe.save`.

## Role

Build adaptive ISDF into `./SAVE/` for the supercell case `Si8_sc` (`isdf_source_dir`).

## Run

```matlab
cd example
run_sc          % UC then SC
% or UC only:
run_sc('uc')
```

Groundstate lives in `./qe.save/` (~1.3 MB).
