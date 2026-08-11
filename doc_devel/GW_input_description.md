<style>
a { color: #0077cc; text-decoration: none; }
a:hover { text-decoration: underline; }
</style>

# GWOptions Input File Description

This document describes all supported input blocks and parameters used to configure the GWOptions framework.
Parameter names and defaults follow `util/allowed_param_list.m` and `util/default_param_values.m`.

**Default energy unit is Ry, ALWAYS.**
**This could be really delicate!!!! Since codes use different unit for their cutoffs.**

ISDF type suffixes:

| Suffix | Pair | Meaning |
|--------|------|---------|
| `type1` | `vc` | valence–conduction |
| `type2` | `vn` | valence–valence (occupied products for exchange-like kernels) |
| `type3` | `nn` | conduction–conduction |

---

## Alphabetical Parameter Index

### &CONTROL
<a href="#isgw">isgw</a> |
<a href="#isbse">isbse</a> |
<a href="#enablek">enable_k_points</a> |
<a href="#output_dir">output_dir</a> |
<a href="#prefix">prefix</a> |
<a href="#groundstate_dir">groundstate_dir</a> |
<a href="#groundstate_type">groundstate_type</a> |
<a href="#storage_dir">storage_dir</a> |
<a href="#outfile">outfile</a> |
<a href="#log_level">log_level</a> |
<a href="#log_file">log_file</a>

### &SYSTEM
<a href="#number_bands_in_summation">number_bands_in_summation</a> |
<a href="#energy_band_index_min">energy_band_index_min</a> |
<a href="#energy_band_index_max">energy_band_index_max</a>

### &CUTOFFS
<a href="#coulomb_truncation_method">coulomb_truncation_method</a> |
<a href="#coulomb_truncation_parameter">coulomb_truncation_parameter</a> |
<a href="#coulomb_cutoff">coulomb_cutoff</a> |
<a href="#density_cutoff">density_cutoff</a>

### &FREQUENCY
<a href="#frequency_dependence">frequency_dependence</a> |
<a href="#frequency_dependence_method">frequency_dependence_method</a> |
<a href="#frequency_low_cutoff">frequency_low_cutoff</a> |
<a href="#broadening">broadening</a> |
<a href="#delta_frequency">delta_frequency</a> |
<a href="#number_imaginary_freqs">number_imaginary_freqs</a> |
<a href="#eta">eta</a> |
<a href="#cd_integration_method">cd_integration_method</a> |
<a href="#cd_integration_parameter">cd_integration_parameter</a> |
<a href="#cd_residual_method">cd_residual_method</a>

### &ISDF
<a href="#isisdf">isisdf</a> |
<a href="#compute_vc">compute_vc</a> |
<a href="#compute_vn">compute_vn</a> |
<a href="#compute_nn">compute_nn</a> |
<a href="#isdf_ratio_type1">isdf_ratio_type1</a> |
<a href="#isdf_ratio_type2">isdf_ratio_type2</a> |
<a href="#isdf_ratio_type3">isdf_ratio_type3</a> |
<a href="#adaptive_threshold_type1">adaptive_threshold_type1</a> |
<a href="#adaptive_threshold_type2">adaptive_threshold_type2</a> |
<a href="#adaptive_threshold_type3">adaptive_threshold_type3</a> |
<a href="#adaptive_num_add_type1">adaptive_num_add_type1</a> |
<a href="#adaptive_num_add_type2">adaptive_num_add_type2</a> |
<a href="#adaptive_num_add_type3">adaptive_num_add_type3</a> |
<a href="#adaptive_candidate_ratio_type1">adaptive_candidate_ratio_type1</a> |
<a href="#adaptive_candidate_ratio_type2">adaptive_candidate_ratio_type2</a> |
<a href="#adaptive_candidate_ratio_type3">adaptive_candidate_ratio_type3</a> |
<a href="#adaptive_max_add_frac_type1">adaptive_max_add_frac_type1</a> |
<a href="#adaptive_max_add_frac_type2">adaptive_max_add_frac_type2</a> |
<a href="#adaptive_max_add_frac_type3">adaptive_max_add_frac_type3</a> |
<a href="#adaptive_max_cond_type1">adaptive_max_cond_type1</a> |
<a href="#adaptive_max_cond_type2">adaptive_max_cond_type2</a> |
<a href="#adaptive_max_cond_type3">adaptive_max_cond_type3</a> |
<a href="#adaptive_use_cond_guard">adaptive_use_cond_guard</a> |
<a href="#adaptive_batch_size">adaptive_batch_size</a> |
<a href="#exxmethod_type1">exxmethod_type1</a> |
<a href="#exxmethod_type2">exxmethod_type2</a> |
<a href="#exxmethod_type3">exxmethod_type3</a> |
<a href="#exxmethod">exxmethod</a> |
<a href="#seed">seed</a> |
<a href="#init">init</a> |
<a href="#weight">weight</a> |
<a href="#sys">sys</a> |
<a href="#is_helper">is_helper</a> |
<a href="#validate_hf">validate_hf</a> |
<a href="#inv_strategy">inv_strategy</a> |
<a href="#inv_param">inv_param</a> |
<a href="#inv_ratio">inv_ratio</a> |
<a href="#auto_inv_param">auto_inv_param</a> |
<a href="#order">order</a>

### &COHSEX
<a href="#exact_ch">exact_ch</a> |
<a href="#ex_use_which_isdf">ex_use_which_isdf</a>

### &SUPERCELL
<a href="#is_supercell">is_supercell</a> |
<a href="#use_sc_isdf">use_sc_isdf</a> |
<a href="#sc_adaptive">sc_adaptive</a> |
<a href="#k1">k1</a> |
<a href="#k2">k2</a> |
<a href="#k3">k3</a> |
<a href="#isdf_source_dir">isdf_source_dir</a>

### &FORMAL
<a href="#wf_seed">wf_seed</a> |
<a href="#eo_e0">eo_e0</a> |
<a href="#eo_delta">eo_delta</a> |
<a href="#run_qp">run_qp</a> |
<a href="#enable_profile">enable_profile</a>

---

## Namelist

### Namelist: &CONTROL

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="isgw"></a>`isgw` | int | No | `1` | Enable GW calculation |
| <a name="isbse"></a>`isbse` | int | No | `0` | Enable BSE calculation |
| <a name="enablek"></a>`enable_k_points` | int / bool | No | `0` | Enable multi-k calculation |
| <a name="output_dir"></a>`output_dir` | string | No | `'./'` | Output directory |
| <a name="prefix"></a>`prefix` | string | No | `'QP'` | File name prefix |
| <a name="groundstate_dir"></a>`groundstate_dir` | string | Yes | *none* | Path to ground-state data directory |
| <a name="groundstate_type"></a><a href="#appendix-sys-freq">`groundstate_type`</a> | string | Yes | `'kssolv'` | Ground-state source (`kssolv`, `qe`, `formal`, …) |
| <a name="storage_dir"></a>`storage_dir` | string | No | `'./QP.save/'` | Directory for intermediate quantities |
| <a name="outfile"></a>`outfile` | string | No | `'GWoutput'` | Output energy file name |
| <a name="log_level"></a>`log_level` | int | No | `1` | Logging verbosity |
| <a name="log_file"></a>`log_file` | string | No | `[]` | Log file path; empty → console |

---

### Namelist: &SYSTEM

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="number_bands_in_summation"></a>`number_bands_in_summation` | int | <a href="#appendix-sys-freq">system based</a> | `-1` | Band count used in summation (`-1` → fill from ground state) |
| <a name="energy_band_index_min"></a>`energy_band_index_min` | int | <a href="#appendix-sys-freq">system based</a> | `-1` | Minimum band index for QP energies |
| <a name="energy_band_index_max"></a>`energy_band_index_max` | int | <a href="#appendix-sys-freq">system based</a> | `-1` | Maximum band index for QP energies |

---

### Namelist: &CUTOFFS

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="coulomb_truncation_method"></a><a href="#appendix-trunc">`coulomb_truncation_method`</a> | int | No | `2` | Coulomb truncation scheme |
| <a name="coulomb_truncation_parameter"></a>`coulomb_truncation_parameter` | float | No | `5.0` | Truncation parameter (Ry; e.g. sphere radius for method `2`) |
| <a name="coulomb_cutoff"></a>`coulomb_cutoff` | float | No | `-1.0` | Coulomb / G-vector cutoff (Ry). `-1` → <a href="#appendix-sys-freq">from ground state</a>. **Sensitive; treat carefully.** |
| <a name="density_cutoff"></a>`density_cutoff` | float | No | `-1` | Density cutoff (Ry). `-1` → <a href="#appendix-sys-freq">from ground state</a> |

---

### Namelist: &FREQUENCY

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="frequency_dependence"></a>`frequency_dependence` | int | No | `0` | Frequency treatment (`0` = COHSEX / static-like; nonzero selects frequency-dependent paths) |
| <a name="frequency_dependence_method"></a>`frequency_dependence_method` | int | No | `2` | Method for frequency dependence (used when frequency-dependent) |
| <a name="frequency_low_cutoff"></a>`frequency_low_cutoff` | float | No | `-1.0` | Low-frequency cutoff (Ry). `-1` → <a href="#appendix-sys-freq">from system</a> |
| <a name="broadening"></a>`broadening` | float | No | `0.018376` | Broadening for frequency-dependent paths (eV). Callers (e.g. `gw.fullfreq_cd_core_Gamma`) convert to Ry before `isdf.gen_Kq_Gamma` / dense chi |
| <a name="delta_frequency"></a>`delta_frequency` | float | No | `0.146997295433511` | Frequency grid step (Ry) |
| <a name="number_imaginary_freqs"></a>`number_imaginary_freqs` | int | No | `15` | Number of imaginary frequencies |
| <a name="eta"></a>`eta` | float | No | `1e-4` | Broadening |
| <a name="cd_integration_method"></a>`cd_integration_method` | int | No | `0` | Contour-deformation integration method |
| <a name="cd_integration_parameter"></a>`cd_integration_parameter` | float | No | `2.0` | CD integration parameter (Ry) |
| <a name="cd_residual_method"></a>`cd_residual_method` | int | No | `0` | CD residual / resolution method |

---

### Namelist: &ISDF

Target rank scale: \(N_\mu \approx \texttt{isdf\_ratio\_type*} \sqrt{N_1 N_2}\) (adaptive paths also fold in \(\sqrt{n_{\mathrm{bz}}}\)).

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="isisdf"></a>`isisdf` | int / bool | No | `1` | Enable ISDF |
| <a name="compute_vc"></a>`compute_vc` | bool | No | `true` | Build / adapt ISDF for type1 (`vc`) |
| <a name="compute_vn"></a>`compute_vn` | bool | No | `true` | Build / adapt ISDF for type2 (`vn`) |
| <a name="compute_nn"></a>`compute_nn` | bool | No | `true` | Build / adapt ISDF for type3 (`nn`) |
| <a name="isdf_ratio_type1"></a>`isdf_ratio_type1` | float | No | `8.0` | ISDF ratio for `vc` (must be &gt; 0) |
| <a name="isdf_ratio_type2"></a>`isdf_ratio_type2` | float | No | `8.0` | ISDF ratio for `vn` (must be &gt; 0) |
| <a name="isdf_ratio_type3"></a>`isdf_ratio_type3` | float | No | `8.0` | ISDF ratio for `nn` (must be &gt; 0) |
| <a name="adaptive_threshold_type1"></a>`adaptive_threshold_type1` | float | No | `2e-4` | Adaptive loss threshold for `vc`. Also used to pick single vs double adaptive backend (see note below) |
| <a name="adaptive_threshold_type2"></a>`adaptive_threshold_type2` | float | No | `2e-6` | Adaptive loss threshold for `vn` |
| <a name="adaptive_threshold_type3"></a>`adaptive_threshold_type3` | float | No | `2e-6` | Adaptive loss threshold for `nn` |
| <a name="adaptive_num_add_type1"></a>`adaptive_num_add_type1` | int | No | `16` | Points added per adaptive iteration (`vc`) |
| <a name="adaptive_num_add_type2"></a>`adaptive_num_add_type2` | int | No | `16` | Points added per adaptive iteration (`vn`) |
| <a name="adaptive_num_add_type3"></a>`adaptive_num_add_type3` | int | No | `16` | Points added per adaptive iteration (`nn`) |
| <a name="adaptive_candidate_ratio_type1"></a>`adaptive_candidate_ratio_type1` | float | No | `2.0` | Candidate-pool size factor vs `num_add` (`vc`) |
| <a name="adaptive_candidate_ratio_type2"></a>`adaptive_candidate_ratio_type2` | float | No | `2.0` | Candidate-pool size factor vs `num_add` (`vn`) |
| <a name="adaptive_candidate_ratio_type3"></a>`adaptive_candidate_ratio_type3` | float | No | `2.0` | Candidate-pool size factor vs `num_add` (`nn`) |
| <a name="adaptive_max_add_frac_type1"></a>`adaptive_max_add_frac_type1` | float | No | `1.0` | Reserved: intended cap factor on adaptive \(N_\mu\) for `vc` (accepted in input; not currently applied in adaptive code) |
| <a name="adaptive_max_add_frac_type2"></a>`adaptive_max_add_frac_type2` | float | No | `1.0` | Reserved: same for `vn` |
| <a name="adaptive_max_add_frac_type3"></a>`adaptive_max_add_frac_type3` | float | No | `1.0` | Reserved: same for `nn` |
| <a name="adaptive_max_cond_type1"></a>`adaptive_max_cond_type1` | float | No | `1e12` | Max allowed condition number in adaptive Schur updates (`vc`) |
| <a name="adaptive_max_cond_type2"></a>`adaptive_max_cond_type2` | float | No | `1e12` | Same for `vn` |
| <a name="adaptive_max_cond_type3"></a>`adaptive_max_cond_type3` | float | No | `1e12` | Same for `nn` |
| <a name="adaptive_use_cond_guard"></a>`adaptive_use_cond_guard` | bool | No | `true` | Enable condition-number guard during adaptive ISDF |
| <a name="adaptive_batch_size"></a>`adaptive_batch_size` | int | No | `64` | Batch size for adaptive weight / residual evaluation |
| <a name="exxmethod_type1"></a>`exxmethod_type1` | string | No | `''` → `'pseudo'` | Initial interpolation-point method for `vc`. Empty filled from `exxmethod`, else `'pseudo'` |
| <a name="exxmethod_type2"></a>`exxmethod_type2` | string | No | `''` → `'pseudo'` | Same for `vn` |
| <a name="exxmethod_type3"></a>`exxmethod_type3` | string | No | `''` → `'pseudo'` | Same for `nn` |
| <a name="exxmethod"></a>`exxmethod` | string | No | `''` → `'pseudo'` | Global fallback for empty `exxmethod_type*` (e.g. `'pseudo'`, `'coarse'`, `'qrcp'`, `'kmeans'`). Non-default methods emit a warning in `gen_coeff` |
| <a name="seed"></a>`seed` | int | No | `0` | RNG seed (legacy / helper paths) |
| <a name="init"></a>`init` | string | No | `'random'` | Initialization strategy (legacy / helper paths) |
| <a name="weight"></a>`weight` | string | No | `'add'` | Weighting strategy label (legacy; not the adaptive residual weight) |
| <a name="sys"></a>`sys` | struct | No | `[]` | Optional system handle; filled from loaded ground state when empty |
| <a name="is_helper"></a>`is_helper` | bool | No | `false` | Enable helper quantities / side paths |
| <a name="validate_hf"></a>`validate_hf` | bool | No | `false` | After ISDF build, run HF-style validation reports |
| <a name="inv_strategy"></a>`inv_strategy` | string | No | `'dir'` | Strategy label for Coulomb / \(\tilde V_q\) inversion bookkeeping |
| <a name="inv_param"></a>`inv_param` | float | No | `1e-12` | SVD / singular-value cutoff (`s_cut`) for \(\tilde V_q\) construction |
| <a name="inv_ratio"></a>`inv_ratio` | float | No | `0.75` | Exponent split for truncated singular values in \(\tilde V_q\) / \(K_q\) (clamped to \([0,1]\)) |
| <a name="auto_inv_param"></a>`auto_inv_param` | bool | No | `false` | Auto-tune `inv_param` when supported |
| <a name="order"></a>`order` | int | No | `1` | Expansion / ordering control for \(\tilde V_q\)-related paths |

**Adaptive backend note.** Routing uses `constant_map().ADAPTIVE_THRESHOLD_CUTOFF` (`1e-5`): `vc` uses `adaptive_single`; `vn` / `nn` use `adaptive_double` when the corresponding threshold is \(\le\) that cutoff, otherwise `adaptive_single`.

---

### Namelist: &COHSEX

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="exact_ch"></a>`exact_ch` | bool | No | `false` | Use exact COH (Coulomb-hole) path when available. Forced `true` for `groundstate_type='formal'` |
| <a name="ex_use_which_isdf"></a>`ex_use_which_isdf` | string | No | `'vn'` | Which ISDF slot feeds exchange: `'vn'` or `'nn'` |

---

### Namelist: &SUPERCELL

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="is_supercell"></a>`is_supercell` | bool | No | `false` | Mark run as supercell workflow |
| <a name="use_sc_isdf"></a>`use_sc_isdf` | bool | No | `false` | Use supercell ISDF import / SC_ISDF path. Requires `is_supercell` and a non-empty `isdf_source_dir` |
| <a name="sc_adaptive"></a>`sc_adaptive` | bool | No | `false` | Run adaptive refinement inside SC_ISDF. Requires `use_sc_isdf=true` |
| <a name="k1"></a>`k1` | int | No | `1` | Supercell replication along \(a_1^*\) |
| <a name="k2"></a>`k2` | int | No | `1` | Supercell replication along \(a_2^*\) |
| <a name="k3"></a>`k3` | int | No | `1` | Supercell replication along \(a_3^*\) |
| <a name="isdf_source_dir"></a>`isdf_source_dir` | string | No\* | `''` | Directory of source ISDF / checkpoint data (\*required when `use_sc_isdf=true`) |

---

### Namelist: &FORMAL

Used when `groundstate_type='formal'` (synthetic / benchmark ground state).

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| <a name="wf_seed"></a>`wf_seed` | int | No | `42` | RNG seed for synthetic wavefunctions |
| <a name="eo_e0"></a>`eo_e0` | float | No | `-2.0` | Base orbital energy (Ry): \(E_o(i_b)=e_0+(i_b-1)\Delta\) |
| <a name="eo_delta"></a>`eo_delta` | float | No | `0.05` | Orbital energy spacing (Ry) |
| <a name="run_qp"></a>`run_qp` | bool | No | `true` | Run QP / COHSEX package after setup |
| <a name="enable_profile"></a>`enable_profile` | bool | No | `false` | Write MATLAB profiler HTML under the case directory |

---

## 📎 Appendix

### <a name="appendix-sys-freq"></a> System-dependent parameters

The default value of `frequency_low_cutoff`, `coulomb_cutoff`, `density_cutoff` are not static when left at `-1`.
They depend on the ground-state information in `groundstate_dir`.

- `coulomb_cutoff` defaults to the wavefunction cutoff of the ground state; `density_cutoff` defaults to twice that cutoff.
- `frequency_low_cutoff` depends on the band energies:
  - Insulators: larger of (HOMO − lowest occupied) and (highest unoccupied computed − LUMO).
  - Metals: highest unoccupied computed − lowest occupied.

> *TODO: Insert algorithm description for estimating frequency_low_cutoff based on sys structure.*

### <a name="appendix-trunc"></a> Coulomb truncation

Value of `coulomb_truncation_method` and corresponding method:

- `0` — no truncation (3D)
- `2` — 0D spherical truncation
- `4` — 1D cell wire truncation **Ongoing**
- `5` — 0D cell box truncation **Ongoing**
- `6` — 2D slab truncation **Ongoing**
- `7` — supercell truncation (3D), experimental **Ongoing**

---

### <a name="gs_input"></a> Requirements for different groundstate types

#### groundstate_type = `kssolv`

**User should call a helper function `save_groundstate_to_GWformat` at the end of your KSSOLV simulation to save these in a standard format. Users should copy `save_groundstate_to_GWformat.m` to their working directory and change their demo in the following manner.**

**A sample demo**
```matlab
% Prepare KSSOLV-scf input
[mol,H,X0,info] = scf(mol, options_scf);
output_dir = './';
save_groundstate_to_GWformat(mol, H, X0, info, output_dir);
```

#### groundstate_type = `qe`

> TODO: complete this part

#### groundstate_type = `formal`

Requires `&SUPERCELL use_sc_isdf = .true.` (and a valid `isdf_source_dir`). Formal mode also forces `validate_hf = false` and `exact_ch = true` in post-processing of defaults.

---

## Notes

- All parameter names are case-insensitive.
- Input follows a namelist-style format (`&BLOCK ... END &BLOCK` or `/`).
- Comments starting with `!` or `#` are ignored.
- Source of truth for allowed keys / defaults: `input/allowed_param_list.m`, `input/default_param_values.m`.
