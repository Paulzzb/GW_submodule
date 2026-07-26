# test_SC Validation Report — Si (vc + vn)

Generated: 2026-06-11 16:24:31

## Configuration

- `compute_vc = .true.` and `compute_vn = .true.` in Si2 / Si8 / Si16 `test` files
- `energy_band_index_max` = 2 × occupied bands (Si2: 8, Si8: 32, Si16: 64)
- LiH_666 stopped / not included

## HF Validation Summary

| Case | Type | id | nisdf | ‖ΔE‖₂/‖E_HF‖₂ |
|------|------|-----|-------|----------------|
| Si2 | vc | 2 | 19 | 6.42760221e-07 |
| Si2 | vn | 4 | 26 | 2.16195033e-06 |
| Si8_from_Si2 | vc | 1 | 76 | 5.43428518e-01 |
| Si8_from_Si2 | vn | 2 | 104 | 1.00505878e+00 |
| Si16_from_Si2 | vc | 1 | 152 | 5.67230507e-01 |
| Si16_from_Si2 | vn | 2 | 208 | 3.97002021e-01 |

## Observations

- **Si2 unit cell:** both `vc` and `vn` adaptive ISDF reach machine-precision HF agreement (relative L2 norm ≲ 10⁻⁶).
- **Si8 / Si16 SC_ISDF:** `vc` SC replication shows smaller relative L2 error than `vn` on the same supercell grid in this run; both remain far from unit-cell accuracy.
- **Si16** `numerical_cond_report` includes both `vc_sc` and `vn_sc` `gen_tildeVq` blocks; Si2/Si8 reports on disk may only reflect partial vn entries (see appendices).

---

# Appendices


## Appendix A: Si2 (unit cell, adaptive)

### A.1 `Si2/numerical_cond_report.txt`

```
=== isdf numerical_cond_report ===
Generated: 2026-06-11 16:20:04
pwd: /data-storage/home/huwei/zzb/GW_double_test/test_profile/test_SC/Si2

--- gen_tildeVq ---
id: 2
desc: vc
s_cut (inv_param): 1.0000e-12

  iqibz: 1
  CCHq cond_number: 1.6010e+04
  CCHq sigma_max: 3.7456e-03
  CCHq sigma_min: 2.3396e-07
  CCHq fro_norm: 6.3539e-03
  CCHq l_keep (SVD truncate): 16
  CCHq sum_truncated_sigma: 6.2025e-19
  CCHq s_cut: 1.0000e-12
  MCHq fro_norm: 5.2975e-02

--- gen_tildeVq ---
id: 4
desc: vn
s_cut (inv_param): 1.0000e-12

  iqibz: 1
  CCHq cond_number: 4.4058e+05
  CCHq sigma_max: 1.4830e-02
  CCHq sigma_min: 3.3660e-08
  CCHq fro_norm: 2.6844e-02
  CCHq l_keep (SVD truncate): 26
  CCHq sum_truncated_sigma: 0.0000e+00
  CCHq s_cut: 1.0000e-12
  MCHq fro_norm: 1.6792e-01
```

### A.2 `Si2/isdf_validate_HF_id2.txt (vc)`

```
=== ISDF validate HF report (isdf id = 2) ===
Generated: 2026-06-11 16:20:36
description: vc
interp_scheme: adaptive
nisdf: 19

--- Band-resolved summary (E_HF vs E_HF_ISDF) ---
Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.

<strong>Item</strong>                                           <strong>Result</strong>                  
    <strong>_______________________________________________</strong>    <strong>_________________________________________</strong>

    "sum(E_HF)"                                        "6.21179023e-01"                         
    "sum(E_HF_ISDF)"                                   "6.21179023e-01"                         
    "sum(E_HF) - sum(E_HF_ISDF)"                       "-2.83953971e-11"                        
    "mean(E_HF) (over nb*nibz*nspin cells)"            "7.76473779e-03"                         
    "mean(E_HF_ISDF) (over nb*nibz*nspin cells)"       "7.76473779e-03"                         
    "sum(|E_HF - E_HF_ISDF|)"                          "3.13942544e-07"                         
    "max(|E_HF - E_HF_ISDF|)"                          "1.56957074e-07"                         
    "||vec(E_HF - E_HF_ISDF)||_2"                      "2.06640979e-07"                         
    "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2"    "6.42760221e-07"                         
    "E_HF minimum"                                     "0.00000000e+00 at (ib=1, ik=1, ispin=1)"
    "E_HF maximum"                                     "2.27179931e-01 at (ib=8, ik=1, ispin=1)"
    "E_HF_ISDF minimum"                                "0.00000000e+00 at (ib=1, ik=1, ispin=1)"
    "E_HF_ISDF maximum"                                "2.27179931e-01 at (ib=8, ik=1, ispin=1)"

--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---
preview_rows budget: 1

<strong>Ex_t</strong>      <strong>Ex_ISDF</strong>        <strong>Diff</strong>         <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
    <strong>________</strong>    <strong>________</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    0.043633    0.043633    -4.8715e-09    4.8715e-09       -1.1165e-07             1.1165e-07

--- Statistics (Mean / Std / Var / Max), ob layer ---
N (total ob samples) = 16
N_rel (samples with |Ex_t| >= tol for relative stats) = 16

<strong>Ex_t</strong>        <strong>Ex_ISDF</strong>         <strong>Diff</strong>         <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
            <strong>__________</strong>    <strong>__________</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    <strong>Mean</strong>      0.038824      0.038824    -1.7747e-12    3.2834e-08       -4.7448e-08             9.1297e-07      
    <strong>Std </strong>      0.014122      0.014122     4.8998e-08    3.5367e-08        1.3325e-06             9.4286e-07      
    <strong>Var </strong>    0.00019944    0.00019944     2.4008e-15    1.2508e-15        1.7757e-12             8.8898e-13      
    <strong>Max </strong>      0.064696      0.064696     1.1102e-07    1.1102e-07        2.5443e-06             2.7053e-06

Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).

--- Global energy sums (accumulated over ob paths) ---
Sum Ex_t^2 (direct):      2.710803953991e-02
Sum Ex_ISDF^2:            2.710803829197e-02
Sum (Ex_t - Ex_ISDF)^2:   3.601158667794e-14

=== end of report ===
```

### A.3 `Si2/isdf_validate_HF_id4.txt (vn)`

```
=== ISDF validate HF report (isdf id = 4) ===
Generated: 2026-06-11 16:20:37
description: vn
interp_scheme: adaptive
nisdf: 26

--- Band-resolved summary (E_HF vs E_HF_ISDF) ---
Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.

<strong>Item</strong>                                           <strong>Result</strong>                  
    <strong>_______________________________________________</strong>    <strong>_________________________________________</strong>

    "sum(E_HF)"                                        "6.13756826e+00"                         
    "sum(E_HF_ISDF)"                                   "6.13756826e+00"                         
    "sum(E_HF) - sum(E_HF_ISDF)"                       "-3.43492346e-10"                        
    "mean(E_HF) (over nb*nibz*nspin cells)"            "7.67196032e-02"                         
    "mean(E_HF_ISDF) (over nb*nibz*nspin cells)"       "7.67196032e-02"                         
    "sum(|E_HF - E_HF_ISDF|)"                          "9.52391650e-06"                         
    "max(|E_HF - E_HF_ISDF|)"                          "4.50142605e-06"                         
    "||vec(E_HF - E_HF_ISDF)||_2"                      "6.00606551e-06"                         
    "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2"    "2.16195033e-06"                         
    "E_HF minimum"                                     "0.00000000e+00 at (ib=9, ik=1, ispin=1)"
    "E_HF maximum"                                     "1.40276350e+00 at (ib=4, ik=1, ispin=1)"
    "E_HF_ISDF minimum"                                "0.00000000e+00 at (ib=9, ik=1, ispin=1)"
    "E_HF_ISDF maximum"                                "1.40276800e+00 at (ib=4, ik=1, ispin=1)"

--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---
preview_rows budget: 1

<strong>Ex_t</strong>     <strong>Ex_ISDF</strong>       <strong>Diff</strong>        <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
    <strong>______</strong>    <strong>_______</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    1.1994    1.1994     2.5363e-10    2.5363e-10        2.1147e-10             2.1147e-10

--- Statistics (Mean / Std / Var / Max), ob layer ---
N (total ob samples) = 32
N_rel (samples with |Ex_t| >= tol for relative stats) = 32

<strong>Ex_t</strong>      <strong>Ex_ISDF</strong>       <strong>Diff</strong>         <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
            <strong>_______</strong>    <strong>_______</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    <strong>Mean</strong>     0.1918     0.1918    -1.0734e-11    3.2169e-07        1.1482e-08             1.4089e-06      
    <strong>Std </strong>    0.40253    0.40253     1.0786e-06    1.0279e-06        2.2289e-06             1.7085e-06      
    <strong>Var </strong>    0.16203    0.16203     1.1633e-12    1.0565e-12         4.968e-12             2.9191e-12      
    <strong>Max </strong>     1.2551     1.2551     3.9553e-06    4.4659e-06        4.3519e-06             6.3083e-06

Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).

--- Global energy sums (accumulated over ob paths) ---
Sum Ex_t^2 (direct):      6.200208544926e+00
Sum Ex_ISDF^2:            6.200208354828e+00
Sum (Ex_t - Ex_ISDF)^2:   3.606280313880e-11

=== end of report ===
```


## Appendix B: Si8_from_Si2 (SC 2×2×1)

### B.1 `Si8_from_Si2/numerical_cond_report.txt`

```
=== isdf numerical_cond_report ===
Generated: 2026-06-11 16:21:04
pwd: /data-storage/home/huwei/zzb/GW_double_test/test_profile/test_SC/Si8_from_Si2

--- gen_tildeVq ---
id: 1
desc: vc_sc_2_2_1
s_cut (inv_param): 1.0000e-12

  iqibz: 1
  CCHq cond_number: 3.7317e+03
  CCHq sigma_max: 6.8705e-04
  CCHq sigma_min: 1.8411e-07
  CCHq fro_norm: 1.3092e-03
  CCHq l_keep (SVD truncate): 16
  CCHq sum_truncated_sigma: 1.7385e-18
  CCHq s_cut: 1.0000e-12
  MCHq fro_norm: 1.3286e-02

--- gen_tildeVq ---
id: 2
desc: vn_sc_2_2_1
s_cut (inv_param): 1.0000e-12

  iqibz: 1
  CCHq cond_number: 1.2969e+09
  CCHq sigma_max: 2.3226e-03
  CCHq sigma_min: 1.7909e-12
  CCHq fro_norm: 4.7429e-03
  CCHq l_keep (SVD truncate): 32
  CCHq sum_truncated_sigma: 6.5149e-18
  CCHq s_cut: 1.0000e-12
  MCHq fro_norm: 3.7398e-02
```

### B.2 `Si8_from_Si2/isdf_validate_HF_id1.txt (vc)`

```
=== ISDF validate HF report (isdf id = 1) ===
Generated: 2026-06-11 16:21:39
description: vc_sc_2_2_1
interp_scheme: adaptive_sc
nisdf: 76

--- Band-resolved summary (E_HF vs E_HF_ISDF) ---
Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.

<strong>Item</strong>                                            <strong>Result</strong>                  
    <strong>_______________________________________________</strong>    <strong>__________________________________________</strong>

    "sum(E_HF)"                                        "2.25868508e+01"                          
    "sum(E_HF_ISDF)"                                   "1.38991977e+01"                          
    "sum(E_HF) - sum(E_HF_ISDF)"                       "8.68765305e+00"                          
    "mean(E_HF) (over nb*nibz*nspin cells)"            "2.82335635e-01"                          
    "mean(E_HF_ISDF) (over nb*nibz*nspin cells)"       "1.73739971e-01"                          
    "sum(|E_HF - E_HF_ISDF|)"                          "1.09583539e+01"                          
    "max(|E_HF - E_HF_ISDF|)"                          "8.30018601e-01"                          
    "||vec(E_HF - E_HF_ISDF)||_2"                      "2.51998594e+00"                          
    "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2"    "5.43428518e-01"                          
    "E_HF minimum"                                     "0.00000000e+00 at (ib=33, ik=1, ispin=1)"
    "E_HF maximum"                                     "1.38505995e+00 at (ib=1, ik=1, ispin=1)" 
    "E_HF_ISDF minimum"                                "0.00000000e+00 at (ib=33, ik=1, ispin=1)"
    "E_HF_ISDF maximum"                                "9.54831887e-01 at (ib=5, ik=1, ispin=1)"

--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---
preview_rows budget: 1

<strong>Ex_t</strong>       <strong>Ex_ISDF</strong>      <strong>Diff</strong>      <strong>AbsDiff</strong>    <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
    <strong>_______</strong>    <strong>_________</strong>    <strong>_______</strong>    <strong>_______</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    0.30017    0.0029449    0.29723    0.29723         0.99019                 0.99019

--- Statistics (Mean / Std / Var / Max), ob layer ---
N (total ob samples) = 512
N_rel (samples with |Ex_t| >= tol for relative stats) = 427

<strong>Ex_t</strong>        <strong>Ex_ISDF</strong>       <strong>Diff</strong>        <strong>AbsDiff</strong>     <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
            <strong>_________</strong>    <strong>_________</strong>    <strong>_________</strong>    <strong>_________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    <strong>Mean</strong>     0.052267     0.031891     0.020376     0.036109         0.11159                 0.72693       
    <strong>Std </strong>     0.076473     0.043826     0.074213     0.067947          1.0502                 0.76529       
    <strong>Var </strong>    0.0058481    0.0019207    0.0055075    0.0046167          1.1029                 0.58568       
    <strong>Max </strong>      0.48285      0.19115      0.40619      0.40619         0.99613                  6.5053

Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).

--- Global energy sums (accumulated over ob paths) ---
Sum Ex_t^2 (direct):      3.658771613548e+00
Sum Ex_ISDF^2:            1.255123244418e+00
Sum (Ex_t - Ex_ISDF)^2:   2.524970075685e+00

=== end of report ===
```

### B.3 `Si8_from_Si2/isdf_validate_HF_id2.txt (vn)`

```
=== ISDF validate HF report (isdf id = 2) ===
Generated: 2026-06-11 16:21:41
description: vn_sc_2_2_1
interp_scheme: adaptive_sc
nisdf: 104

--- Band-resolved summary (E_HF vs E_HF_ISDF) ---
Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.

<strong>Item</strong>                                            <strong>Result</strong>                  
    <strong>_______________________________________________</strong>    <strong>__________________________________________</strong>

    "sum(E_HF)"                                        "2.25868508e+01"                          
    "sum(E_HF_ISDF)"                                   "4.42975753e+01"                          
    "sum(E_HF) - sum(E_HF_ISDF)"                       "-2.17107246e+01"                         
    "mean(E_HF) (over nb*nibz*nspin cells)"            "2.82335635e-01"                          
    "mean(E_HF_ISDF) (over nb*nibz*nspin cells)"       "5.53719691e-01"                          
    "sum(|E_HF - E_HF_ISDF|)"                          "2.17107246e+01"                          
    "max(|E_HF - E_HF_ISDF|)"                          "2.39372306e+00"                          
    "||vec(E_HF - E_HF_ISDF)||_2"                      "4.66065711e+00"                          
    "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2"    "1.00505878e+00"                          
    "E_HF minimum"                                     "0.00000000e+00 at (ib=33, ik=1, ispin=1)"
    "E_HF maximum"                                     "1.38505995e+00 at (ib=1, ik=1, ispin=1)" 
    "E_HF_ISDF minimum"                                "0.00000000e+00 at (ib=33, ik=1, ispin=1)"
    "E_HF_ISDF maximum"                                "3.39381063e+00 at (ib=16, ik=1, ispin=1)"

--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---
preview_rows budget: 1

<strong>Ex_t</strong>      <strong>Ex_ISDF</strong>       <strong>Diff</strong>         <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
    <strong>_______</strong>    <strong>_______</strong>    <strong>___________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    0.30017    0.30017    -1.9429e-15    1.9429e-15       -6.4726e-15             6.4726e-15

--- Statistics (Mean / Std / Var / Max), ob layer ---
N (total ob samples) = 512
N_rel (samples with |Ex_t| >= tol for relative stats) = 427

<strong>Ex_t</strong>       <strong>Ex_ISDF</strong>       <strong>Diff</strong>        <strong>AbsDiff</strong>     <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
            <strong>_________</strong>    <strong>________</strong>    <strong>_________</strong>    <strong>_________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    <strong>Mean</strong>     0.052267    0.096607     -0.04434     0.053105         -2.0029                  2.088        
    <strong>Std </strong>     0.076473     0.11825     0.087146     0.082088          3.1286                 3.0723        
    <strong>Var </strong>    0.0058481    0.013984    0.0075944    0.0067384          9.7879                  9.439        
    <strong>Max </strong>      0.48285      1.0628      0.27638       0.7395         0.89148                 23.972

Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).

--- Global energy sums (accumulated over ob paths) ---
Sum Ex_t^2 (direct):      3.658771613548e+00
Sum Ex_ISDF^2:            1.026048896492e+01
Sum (Ex_t - Ex_ISDF)^2:   4.372248392396e+00

=== end of report ===
```


## Appendix C: Si16_from_Si2 (SC 2×2×2)

### C.1 `Si16_from_Si2/numerical_cond_report.txt`

```
=== isdf numerical_cond_report ===
Generated: 2026-06-11 16:18:19
pwd: /data-storage/home/huwei/zzb/GW_double_test/test_profile/test_SC/Si16_from_Si2

--- gen_tildeVq ---
id: 1
desc: vc_sc_2_2_2
s_cut (inv_param): 1.0000e-12

  iqibz: 1
  CCHq cond_number: 7.3585e+01
  CCHq sigma_max: 5.2135e-04
  CCHq sigma_min: 7.0851e-06
  CCHq fro_norm: 8.4103e-04
  CCHq l_keep (SVD truncate): 16
  CCHq sum_truncated_sigma: 2.0392e-18
  CCHq s_cut: 1.0000e-12
  MCHq fro_norm: 8.2499e-03

--- gen_tildeVq ---
id: 2
desc: vn_sc_2_2_2
s_cut (inv_param): 1.0000e-12

  iqibz: 1
  CCHq cond_number: 7.3532e+09
  CCHq sigma_max: 1.1488e-03
  CCHq sigma_min: 1.5623e-13
  CCHq fro_norm: 2.3899e-03
  CCHq l_keep (SVD truncate): 32
  CCHq sum_truncated_sigma: 7.2379e-18
  CCHq s_cut: 1.0000e-12
  MCHq fro_norm: 1.8667e-02
```

### C.2 `Si16_from_Si2/isdf_validate_HF_id1.txt (vc)`

```
=== ISDF validate HF report (isdf id = 1) ===
Generated: 2026-06-11 16:19:00
description: vc_sc_2_2_2
interp_scheme: adaptive_sc
nisdf: 152

--- Band-resolved summary (E_HF vs E_HF_ISDF) ---
Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.

<strong>Item</strong>                                            <strong>Result</strong>                  
    <strong>_______________________________________________</strong>    <strong>__________________________________________</strong>

    "sum(E_HF)"                                        "4.51839615e+01"                          
    "sum(E_HF_ISDF)"                                   "2.04143991e+01"                          
    "sum(E_HF) - sum(E_HF_ISDF)"                       "2.47695624e+01"                          
    "mean(E_HF) (over nb*nibz*nspin cells)"            "5.64799519e-01"                          
    "mean(E_HF_ISDF) (over nb*nibz*nspin cells)"       "2.55179988e-01"                          
    "sum(|E_HF - E_HF_ISDF|)"                          "2.52404509e+01"                          
    "max(|E_HF - E_HF_ISDF|)"                          "8.13182995e-01"                          
    "||vec(E_HF - E_HF_ISDF)||_2"                      "3.62566046e+00"                          
    "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2"    "5.67230507e-01"                          
    "E_HF minimum"                                     "0.00000000e+00 at (ib=65, ik=1, ispin=1)"
    "E_HF maximum"                                     "1.39694331e+00 at (ib=1, ik=1, ispin=1)" 
    "E_HF_ISDF minimum"                                "0.00000000e+00 at (ib=65, ik=1, ispin=1)"
    "E_HF_ISDF maximum"                                "6.90504012e-01 at (ib=5, ik=1, ispin=1)"

--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---
preview_rows budget: 1

<strong>Ex_t</strong>       <strong>Ex_ISDF</strong>      <strong>Diff</strong>      <strong>AbsDiff</strong>    <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
    <strong>_______</strong>    <strong>_________</strong>    <strong>_______</strong>    <strong>_______</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    0.15012    0.0011571    0.14896    0.14896         0.99229                 0.99229

--- Statistics (Mean / Std / Var / Max), ob layer ---
N (total ob samples) = 2048
N_rel (samples with |Ex_t| >= tol for relative stats) = 1821

<strong>Ex_t</strong>        <strong>Ex_ISDF</strong>         <strong>Diff</strong>       <strong>AbsDiff</strong>     <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
            <strong>__________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    <strong>Mean</strong>      0.024368      0.010849       0.01352    0.014894         0.49644                 0.59166       
    <strong>Std </strong>      0.028178      0.015461      0.024961    0.024166         0.46378                 0.33381       
    <strong>Var </strong>    0.00079402    0.00023906    0.00062307    0.000584         0.21509                 0.11143       
    <strong>Max </strong>       0.23535       0.22052       0.23357     0.23357         0.99479                  5.5468

Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).

--- Global energy sums (accumulated over ob paths) ---
Sum Ex_t^2 (direct):      2.529546226394e+00
Sum Ex_ISDF^2:            6.528126595694e-01
Sum (Ex_t - Ex_ISDF)^2:   1.468242281627e+00

=== end of report ===
```

### C.3 `Si16_from_Si2/isdf_validate_HF_id2.txt (vn)`

```
=== ISDF validate HF report (isdf id = 2) ===
Generated: 2026-06-11 16:19:08
description: vn_sc_2_2_2
interp_scheme: adaptive_sc
nisdf: 208

--- Band-resolved summary (E_HF vs E_HF_ISDF) ---
Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.

<strong>Item</strong>                                            <strong>Result</strong>                  
    <strong>_______________________________________________</strong>    <strong>__________________________________________</strong>

    "sum(E_HF)"                                        "4.51839615e+01"                          
    "sum(E_HF_ISDF)"                                   "5.13237755e+01"                          
    "sum(E_HF) - sum(E_HF_ISDF)"                       "-6.13981401e+00"                         
    "mean(E_HF) (over nb*nibz*nspin cells)"            "5.64799519e-01"                          
    "mean(E_HF_ISDF) (over nb*nibz*nspin cells)"       "6.41547194e-01"                          
    "sum(|E_HF - E_HF_ISDF|)"                          "1.42613857e+01"                          
    "max(|E_HF - E_HF_ISDF|)"                          "1.18035379e+00"                          
    "||vec(E_HF - E_HF_ISDF)||_2"                      "2.53758306e+00"                          
    "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2"    "3.97002021e-01"                          
    "E_HF minimum"                                     "0.00000000e+00 at (ib=65, ik=1, ispin=1)"
    "E_HF maximum"                                     "1.39694331e+00 at (ib=1, ik=1, ispin=1)" 
    "E_HF_ISDF minimum"                                "0.00000000e+00 at (ib=65, ik=1, ispin=1)"
    "E_HF_ISDF maximum"                                "2.14710335e+00 at (ib=30, ik=1, ispin=1)"

--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---
preview_rows budget: 1

<strong>Ex_t</strong>      <strong>Ex_ISDF</strong>       <strong>Diff</strong>        <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
    <strong>_______</strong>    <strong>_______</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    0.15012    0.15012    4.7184e-16    4.7184e-16        3.1432e-15             3.1432e-15

--- Statistics (Mean / Std / Var / Max), ob layer ---
N (total ob samples) = 2048
N_rel (samples with |Ex_t| >= tol for relative stats) = 1821

<strong>Ex_t</strong>        <strong>Ex_ISDF</strong>        <strong>Diff</strong>        <strong>AbsDiff</strong>      <strong>Diff_over_abs_Ex_t</strong>    <strong>AbsDiff_over_abs_Ex_t</strong>
            <strong>__________</strong>    <strong>_________</strong>    <strong>__________</strong>    <strong>__________</strong>    <strong>__________________</strong>    <strong>_____________________</strong>

    <strong>Mean</strong>      0.024368     0.027565    -0.0031964      0.011831        -0.087769                0.57332       
    <strong>Std </strong>      0.028178     0.045459      0.026565      0.023998          0.92673                0.73326       
    <strong>Var </strong>    0.00079402    0.0020665    0.00070571    0.00057588          0.85883                0.53767       
    <strong>Max </strong>       0.23535      0.71584      0.065758       0.49808          0.97243                  9.047

Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).

--- Global energy sums (accumulated over ob paths) ---
Sum Ex_t^2 (direct):      2.529546226394e+00
Sum Ex_ISDF^2:            5.164856852610e+00
Sum (Ex_t - Ex_ISDF)^2:   1.317964590313e+00

=== end of report ===
```
