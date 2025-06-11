<style>
a { color: #0077cc; text-decoration: none; }
a:hover { text-decoration: underline; }
</style>

# GWOptions Input File Description

This document describes all supported input blocks and parameters used to configure the GWOptions framework.
Parameters are grouped by block (namelist-style) and include descriptions, expected types, and default values.

**Default energy unit is Ry, ALWAYS.**
**This could be really delicate!!!! Since codes use different unit for their cutoffs.**

---

## Alphabetical Parameter Index

### &CONTROL
<a href="#isgw">isgw</a> |
<a href="#isbse">isbse</a> |
<a href="#output_dir">output_dir</a> |
<a href="#prefix">prefix</a> |
<a href="#groundstate_dir">groundstate_dir</a> |
<a href="#groundstate_type">groundstate_type</a> |
<a href="#storage_dir">storage_dir</a> |
<a href="#outfile">outfile</a> |

### &SYSTEM
<a href="#number_bands_in_summation">number_bands_in_summation</a> |
<a href="#energy_band_index_min">energy_band_index_min</a> |
<a href="#energy_band_index_max">energy_band_index_max</a> |

### &CUTOFFS
<a href="#coulomb_truncation_method">coulomb_truncation_method</a> |
<a href="#coulomb_truncation_parameter">coulomb_truncation_parameter</a> |
<a href="#coulomb_cutoff">coulomb_cutoff</a> |
<a href="#density_cutoff">density_cutoff</a> |

### &FREQUENCY
<a href="#frequency_dependence">frequency_dependence</a> |
<a href="#frequency_dependence_method">frequency_dependence_method</a> |
<a href="#frequency_low_cutoff">frequency_low_cutoff</a> |
<a href="#delta_frequency">delta_frequency</a> |
<a href="#number_imaginary_freqs">number_imaginary_freqs</a> |
<a href="#eta">eta</a> |
<a href="#cd_int_method">cd_int_method</a> |
<a href="#cd_res_method">cd_res_method</a> |

### &ISDF
<a href="#isisdf">isisdf</a> |
<a href="#isdf_ratio">isdf_ratio</a> |
<a href="#isdf_ratio_type1">isdf_ratio_type1</a> |
<a href="#isdf_ratio_type2">isdf_ratio_type2</a> |
<a href="#isdf_ratio_type3">isdf_ratio_type3</a> |
<a href="#exxmethod">exxmethod</a> |
<a href="#seed">seed</a> |
<a href="#init">init</a> |
<a href="#weight">weight</a> |
<a href="#sys">sys</a> |

---

## Namelist

### Namelist: &CONTROL

| Parameter          | Type   | Required | Default      | Description                                      |
|--------------------|--------|----------|---------------|--------------------------------------------------|
| <a name="isgw"></a>`isgw`                | int    | No       | `1`           | Enable GW calculation                           |
| <a name="isbse"></a>`isbse`              | int    | No       | `0`           | Enable BSE calculation                          |
| <a name="output_dir"></a>`output_dir`    | string | No       | `'./'`        | Output directory                                |
| <a name="prefix"></a>`prefix`            | string | No       | `'QP'`        | File name prefix                                |
| <a name="groundstate_dir"></a>`groundstate_dir` | string | Yes | *none*       | Path to ground state data directory             |
| <a name="groundstate_type"></a><a href="#appendix-sys-freq">groundstate_type</a> | string | Yes | `kssolv`       | Software implementing groundstate calculation             |
| <a name="storage_dir"></a>`storage_dir`  | string | No       | `'./QP.save/'`   | Directory for intermediate quantities           |
| <a name="outfile"></a>`outfile`          | string | No       | `'./GWoutput'`  | Output log file name                            |

---

### Namelist: &SYSTEM

| Parameter                             | Type   | Required | Default    | Description                        |
|----------------------------------------|--------|----------|------------|------------------------------------|
| <a name="number_bands_in_summation"></a>`number_bands_in_summation` | int    | <a href="#appendix-sys-freq">system based</a> | `-1`          | Band count in summation            |
| <a name="energy_band_index_min"></a>`energy_band_index_min` | int    | <a href="#appendix-sys-freq">system based</a> | `-1`          | Minimum band index for QP energies            |
| <a name="energy_band_index_max"></a>`energy_band_index_max` | int    | <a href="#appendix-sys-freq">system based</a> | `-1`          | Maximum band index for QP energies            |

---

### Namelist: &CUTOFFS


| Parameter           | Type   | Required | Default | Description                      |
|---------------------|--------|----------|---------|----------------------------------|
| <a name="coulomb_truncation_method"></a>`coulomb_truncation_method`         | <a href="#appendix-trunc"> `spherical` </a> | No | 2 | Truncation method for Coulomb     |
| <a name="coulomb_truncation_parameter"></a>`coulomb_truncation_parameter`   | float  | No | 5.0  | Truncation parameter for Coulomb potential    |
| <a name="coulomb_cutoff"></a>`coulomb_cutoff`   | float | No | 5.0  | Cutoff for Coulomb potential **Could be danger, code is not stable about it**    |
| <a name="density_cutoff"></a>`density_cutoff`   | float  | No | <a href="#appendix-sys-freq">system based</a>  | Cutoff for density               |

---
---

## Namelist: &FREQUENCY


| Parameter                     | Type   | Required | Default             | Description                                |
|-------------------------------|--------|----------|----------------------|--------------------------------------------|
| <a name="frequency_dependence"></a>`frequency_dependence`         | int    | No | `0`       | Whether to compute frequency-dependence   |
| <a name="frequency_dependence_method"></a>`frequency_dependence_method` | int | No | `0`       | Method to compute frequency-dependence    |
| <a name="frequency_low_cutoff"></a>`frequency_low_cutoff`         | float  | No | <a href="#appendix-sys-freq">system based</a> | Low freq cutoff (computed from system info) |
| <a name="delta_frequency"></a>`delta_frequency`                   | float  | No | `2`       | Frequency step                            |
| <a name="number_imaginary_freqs"></a>`number_imaginary_freqs`     | int    | No | `15`      | Number of imaginary frequencies           |
| <a name="eta"></a>`eta`                                           | float  | No | `1e-4`    | Broadening                                |
| <a name="cd_int_method"></a>`cd_int_method`                       | int    | No | `0`       | CD integration method                     |
| <a name="cd_res_method"></a>`cd_res_method`                       | int    | No | `0`       | CD resolution method                      |

---

### Namelist: &ISDF

| Parameter             | Type   | Required | Default        | Description                         |
|------------------------|--------|----------|----------------|-------------------------------------|
| <a name="isisdf"></a>`isisdf`                   | bool   | No | `true`        | Whether to enable ISDF              |
| <a name="isdf_ratio"></a>`isdf_ratio`           | float  | No | `8.0`         | Global ISDF ratio                   |
| <a name="isdf_ratio_type1"</a>`isdf_ratio_type1` | float | No | = `isdf_ratio` | Override ratio for ISDF type 1      |
| <a name="isdf_ratio_type2"></a>`isdf_ratio_type2` | float | No | = `isdf_ratio` | Override ratio for ISDF type 2      |
| <a name="isdf_ratio_type3"></a>`isdf_ratio_type3` | float | No | = `isdf_ratio` | Override ratio for ISDF type 3      |
| <a name="exxmethod"></a>`exxmethod`             | string | No | `'kmeans'`    | Method for EXX                      |
| <a name="seed"></a>`seed`                       | int    | No | `0`           | Random seed                         |
| <a name="init"></a>`init`                       | string | No | `'wrs'`       | Initialization strategy             |
| <a name="weight"></a>`weight`                   | string | No | `'add'`       | Weighting strategy                  |
| <a name="sys"></a>`sys`                         | struct | No | `[]`          | System placeholder object           |
---

## 📎 Appendix

### <a name="appendix-sys-freq"></a> System-dependent parameters 

The default value of `frequency_low_cutoff`, `coulomb_cutoff`, `density_cutoff` are not static.
They depend on the groundstate information provided by user in `groundstate_dir`.
The algorithm to compute it is as follows:

- `coulomb_cutoff` is the wavefunction cutoff of the groundstate, and the `density_cutoff` is twice the wavefunction cutoff.
- `frequency_low_cutoff` depends on the band energies.
  - For insulating systems, the `frequency_low_cutoff` is set to the larger of the following two values: the energy difference between the HOMO and the lowest occupied band, and the energy difference between the highest unoccupied band calculated and the LUMO.
  - For metallic systems, `frequency_low_cutoff` is set to the highest unoccupied band energy calculated minus the lowest occupied band energy.

> *TODO: Insert algorithm description for estimating frequency_low_cutoff based on sys structure.*

### <a name="#appendix-trunc"> </a> Coulomb truncation
> Todo : introduce different Coulomb truncation schemes and their corresponding parameters.
Value of `coulomb_truncation_method` and corresponding method:

- 0,  no truncation (3D)
- 2,  0D spherical truncation
- 4,  1D cell wire truncation **Ongoing**
- 5,  0D cell box truncation **Ongoing**
- 6,  2D slab truncation **Ongoing**
- 7,  supercell truncation (3D), experimental **Ongoing**

---
### <a name="gs_input"></a> Requirements for different groundstate types

####  groundstate_type = `kssolv`

**User should call a helper function `save_groundstate_to_GWformat` at the end of your KSSOLV simulation to save these in a standard format. Users should copy `save_groundstate_to_GWformat.m` to their working directory and change their demo in the following manner.**

**A sample demo**
```matlab
% Prepare KSSOLV-scf input
[mol,H,X0,info] = scf(mol, options_scf);
output_dir = './';
save_groundstate_to_GWformat(mol, H, X0, info, output_dir);
```

####  groundstate_type = `qe`

> TODO: complete this part
---

## Notes

- All parameter names are case-insensitive.
- Input follows a namelist-style format (`&BLOCK ... /`).
- Comments starting with `!` or `#` are ignored.
