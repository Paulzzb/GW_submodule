# Namelist input (user guide)

Use this when editing the case `test` namelist. Only parameters commonly needed for the showcase examples are listed.  
Default energy unit: **Ry** (except `FREQUENCY.broadening`; see table).

中文版：[`../doc_user_ZH/GW_input_description.md`](../doc_user_ZH/GW_input_description.md)

ISDF type suffixes:

| Suffix | Pair | Meaning |
|--------|------|---------|
| `type1` | `vc` | valence–conduction |
| `type2` | `vn` | valence–valence |
| `type3` | `nn` | conduction–conduction |

Repository examples: `examples/cases/*/test`.

---

## `&CONTROL`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `groundstate_dir` | (required) | Ground-state directory, e.g. `'./qe.save'` |
| `groundstate_type` | `'qe'` | Currently only `'qe'` is supported |
| `storage_dir` | `'./SAVE/'` | Intermediate data; examples usually use `'./SAVE'` |
| `output_dir` | `'./'` | Directory for the report and `qp.dat` |
| `prefix` | `'QP'` | Report prefix → `r-<prefix>.log` |
| `log_level` | `1` | Log verbosity |
| `enable_k_points` | `0` | Multi-k; Gamma demos use `.false.` / `0` |
| `isgw` | `1` | Enable GW |
| `isbse` | `0` | BSE (not supported) |

---

## `&SYSTEM`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `number_bands_in_summation` | `-1` | Bands in summation; `-1` → filled from ground state |
| `energy_band_index_min` | `-1` | Lower band index for QP output |
| `energy_band_index_max` | `-1` | Upper band index for QP output |

---

## `&CUTOFFS`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `coulomb_truncation_method` | `2` | `0` none; `2` spherical truncation (common) |
| `coulomb_truncation_parameter` | `5.0` | Truncation parameter (Ry; sphere radius for method `2`) |
| `coulomb_cutoff` | `-1.0` | Coulomb/G cutoff (Ry); `-1` → use GS wavefunction cutoff |
| `density_cutoff` | `-1` | Density cutoff (Ry); `-1` → ~2× WF cutoff; unused today |

Cutoffs are sensitive; compare against the ground state and examples before changing.

---

## `&FREQUENCY`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `frequency_dependence` | `0` | **`-2`** COHSEX (Gamma); **`2`** full-frequency CD (Gamma). Showcase cases use these two |
| `frequency_dependence_method` | `2` | Method choice for full-frequency |
| `frequency_low_cutoff` | `-1.0` | Low-frequency cutoff (Ry); `-1` → estimated from bands |
| `broadening` | `0.25` | Full-frequency broadening (**eV**) |
| `delta_frequency` | `≈0.147` | Real-frequency grid step (Ry) |
| `number_imaginary_freqs` | `15` | Number of imaginary frequencies |
| `cd_integration_parameter` | `2.0` | CD integration parameter (Ry) |

`qp.launcher` routing:

| `frequency_dependence` | Path |
|------------------------|------|
| `-2` | `gw.x_Gamma` + `gw.cohsex_Gamma` |
| `2` | `gw.fullfreq_cd_res_Gamma` + `gw.fullfreq_cd_int_Gamma` |

---

## `&ISDF`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `isisdf` | `1` | `1` ISDF; `0` dense (comparison cases) |
| `compute_vc` / `compute_vn` / `compute_nn` | `true` | Build the corresponding ISDF type |
| `isdf_ratio_type1/2/3` | `8.0` | ISDF ratio per type (must be &gt; 0) |
| `exxmethod` | → `'pseudo'` | Interpolation points: `'pseudo'` / `'qrcp'` / `'kmeans'` / `'coarse'` … |
| `exxmethod_type1/2/3` | empty → follow `exxmethod` | Per-type override |
| `iscauchy` | `false` | Cauchy-related options |
| `validate_hf` | `false` | Write HF-style validation reports (`isdf_report/`) |
| `is_helper` | `false` | Helper / bypass quantities |

Adaptive refinement (`adaptive_*`) has defaults; follow the examples unless you need fine tuning (`examples/cases/*/test`).

Target rank scale: \(N_\mu \approx \texttt{isdf\_ratio} \sqrt{N_1 N_2}\).

---

## `&COHSEX`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `exact_ch` | `false` | Exact COH sum (only with `isisdf=1`) |
| `ex_use_which_isdf` | `'vn'` | Which ISDF for exchange: `'vn'` or `'nn'` |

This block matters mainly when `frequency_dependence = -2`.

---

## Minimal example (QE + full-frequency ISDF)

```text
&CONTROL
  groundstate_dir = '../qe.save',
  groundstate_type = 'qe',
  storage_dir = './SAVE',
  output_dir = './',
  prefix = 'ff_isdf',
  log_level = 1,
END &CONTROL

&SYSTEM
  number_bands_in_summation = 319,
  energy_band_index_min = 1,
  energy_band_index_max = 32,
END &SYSTEM

&CUTOFFS
  coulomb_truncation_method = 2,
  coulomb_truncation_parameter = 5.0,
  coulomb_cutoff = 15.0,
END &CUTOFFS

&FREQUENCY
  frequency_dependence = 2,
  frequency_dependence_method = 2,
  frequency_low_cutoff = 25.0,
  broadening = 0.01876,
  number_imaginary_freqs = 25,
END &FREQUENCY

&ISDF
  isisdf = 1,
  isdf_ratio_type1 = 12.0,
  isdf_ratio_type2 = 24.0,
  isdf_ratio_type3 = 16.0,
  exxmethod = 'qrcp',
END &ISDF
```

More complete examples: `examples/cases/`.
