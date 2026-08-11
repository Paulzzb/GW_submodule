# Output description

This document summarizes the major output files produced by the current QP calculation framework.

Related docs: [`README.md`](README.md), [`class_reference.md`](class_reference.md).  
File-name constants are defined in `util/filename_map.m`.

中文版：[`../doc_user_ZH/output_description.md`](../doc_user_ZH/output_description.md)

---

## Overview

| Location | Files | Produced by |
|----------|-------|-------------|
| `CONTROL.storage_dir` (e.g. `./SAVE`) | `data.mat`, `config.mat`, `relay_stage.mat`, ISDF checkpoints | `input_driver` / `service_driver` / ISDF |
| `CONTROL.output_dir` (e.g. `./`) | `r-<prefix>.log`, `qp.dat` | `display_input_summary` / `qp.fout` |
| Working case directory | `isdf_report/o-ISDF_*` | ISDF report writers |

Default `storage_dir` is `./SAVE/`; examples usually set `./SAVE` explicitly.

---

## 1. `data.mat`

**Path:** `<storage_dir>/data.mat`  
**Writer:** `input_driver`  
**Variable:** `data`

**Description:**  
Ground-state struct loaded via `interfaces/load_groundstate_info` (typically QE `*.save`).

**Purpose:**  
Persist DFT input so later stages can restart without re-reading QE HDF5.

**Format:** MATLAB `-v7.3` (no compression).

---

## 2. `config.mat`

**Path:** `<storage_dir>/config.mat`  
**Writer:** `input_driver`  
**Variable:** `config`

**Description:**  
Parsed namelist plus defaults (`set_default_param_value`).  
Blocks typically include `CONTROL`, `SYSTEM`, `CUTOFFS`, `FREQUENCY`, `ISDF`, `COHSEX`.

**Important:**  
`input_driver` **always rebuilds** `config` from the namelist.  
Editing CUTOFFS / ISDF / FREQUENCY takes effect on the next `input_driver` call even if the stage cache is kept.

Parameter overview: [`GW_input_description.md`](GW_input_description.md).

---

## 3. `relay_stage.mat`

**Path:** `<storage_dir>/relay_stage.mat`  
**Writer:** `relay.save2db` inside `service_driver`  
**Variable:** `relay_stage`

**Description:**  
Cached service-layer objects (expensive, data-derived):

- restored: `system`, `symmetry`, `FFT`, `wave_functions` (and related)
- **not** staged (rebuilt each run): `lattice`, `coulomb`, ISDF

**Purpose:**  
Allow `qp_driver` / a later MATLAB session to `relay.stage_from_db` + `relay.restore` without reconstructing wavefunctions from scratch.

If `relay_stage.mat` exists but `data.mat` is missing, `input_driver` warns and rebuilds from the ground-state directory.

---

## 4. `qp.dat`

**Path:** `<output_dir>/qp.dat` (falls back to `storage_dir` if `output_dir` unset)  
**Writer:** `qp.fout` (called by `qp.launcher`)  
**Class / struct:** result struct `E` — see [`class_reference.md`](class_reference.md)

**Description:**  
Human-readable quasiparticle energy table.

**Units:** All printed energies are in **electronvolts (eV)**.

**Note:**  
`CONTROL.outfile` (default `'GWoutput'`) is currently **not** used by `qp.fout`; the file name is fixed as `qp.dat`.

### Static / COHSEX-style header

Used when `frequency_dependence` is not `2` (e.g. `-2`):

```plaintext
   n         Emf          Eo           X        SX-X          CH         Sig         Vxc        Eqp0
```

One row per band (and per k in the loop order of `qp.fout`).

### Full-frequency header (`frequency_dependence == 2`)

Two lines per band: real parts on the first line, imaginary parts of SX-X / CH / Sig / Eqp0 on the second:

```plaintext
   n         Emf          Eo           X      Re SX-X       Re CH      Re Sig        Vxc     Re Eqp0
                                     Im SX-X       Im CH      Im Sig                Im Eqp0
  13    3.467416    3.467416  -12.445148    1.211312   -0.172058  -11.405895  -10.521226    2.582748 
                                            0.000059    0.000000    0.000059                0.000059 
```

### Column meaning

| Column | Meaning |
|--------|---------|
| `n` | Band index (`energy_band_index_min` … `max`) |
| `Emf` / `Eo` | Mean-field / KS energy (from `system.get().Eo`) |
| `X` | Exact exchange |
| `SX-X` / `CH` | Screened-exchange correction and Coulomb-hole (or full-freq residual / integral) |
| `Sig` | `X + SX-X + CH` |
| `Vxc` | Exchange-correlation reference |
| `Eqp0` | `Eo + Sig - Vxc` |

---

## 5. In-memory result `E`

Not a file by itself, but the primary MATLAB return value:

```matlab
E = qp.launcher(config);
% or
E = qp_driver('./SAVE');
```

Fields: `Eqp`, `Ex`, `Esx_x`, `Ecoh`, `Eqp0` (empty), `fout` (path to `qp.dat`).  
Details: [`class_reference.md`](class_reference.md).

---

## 6. `r-<prefix>.log`

**Path:** `<output_dir>/r-<prefix>.log`  
**Writer:** `service/+output` (opened by `display_input_summary` during `input_driver`)

**Description:**  
Main run report / log (Yambo-style OF naming).

- Default `prefix = 'QP'` → `r-QP.log`
- Examples: `r-Si_gamma.log`, `r-Si2_uc.log`
- Existing files are rotated: `r-Si_gamma_01.log`, …

Verbosity follows `CONTROL.log_level` (`0`–`2`).  
`CONTROL.log_file` is a legacy exclusive-log option; the live report used by current drivers is `r-<prefix>.log`.

---

## 7. ISDF reports and checkpoints

### 7.1 Text reports

**Directory:** `./isdf_report/` (relative to the case working directory)  
**Names** (`filename_map`):

| Pattern | Role |
|---------|------|
| `o-ISDF_cond` | Conditioning / condition-number report |
| `o-ISDF_HF_id%d` | Hartree–Fock validation report for ISDF id |
| `o-ISDF_adaptive_id%d` | Adaptive ISDF progress / diagnostics |

Enabled when ISDF is on (`ISDF.isisdf`) and the corresponding validation / adaptive paths run.

### 7.2 Adaptive checkpoints

**Path:** `<storage_dir>/isdf_adaptive_checkpoint_<desc>_idN.mat`  
(e.g. `vc` / `vn` / `nn`)

Used to resume expensive adaptive ISDF constructions.

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-11 | ZZ | Split into `doc_user_ZH` / `doc_user_EN` |
| 2026-08-08 | Zhengbang | Rewrite for SAVE / qp.dat / r-*.log layout |
