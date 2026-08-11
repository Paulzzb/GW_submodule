# Introduction and reference of classes / data structures

This document describes the main data carriers used by the QP codebase:

- Ground-state struct `data` (on disk as `data.mat`)
- Service-layer managers (`system`, `lattice`, `coulomb`, …)
- `@gvec` under `src/@gvec` (G-grid construction helper)
- Result struct `E` from `qp.launcher` (with text table `qp.dat`)

Runtime physics kernels read from **service managers**.

For convenience of developers, important fields are summarized below.

---

---

## `data` (ground-state struct)

### A. Introduction

`data` is the plain MATLAB struct produced by `load_groundstate_info` and saved as `SAVE/data.mat`.  
It is the Stage-I / Stage-II carrier of DFT information before service objects are built.

Sources (`CONTROL.groundstate_type`):

- `'qe'` — load QE `*.save` folder via `load_qe_from_folder`
- `'kssolv'` — *(currently unsupported)* load `<dir>/groundstate.mat` (variable `groundstate`)

### B. Properties list

[rhor](#rhor) | [Vxc](#vxc) | [ev](#ev) | [psig](#psig) | [occupation](#occupation) | [reciprocal_grid_info](#reciprocal_grid_info) | [nkibz](#nkibz) | [kibz](#kibz) | [kweight](#kweight) | [nspin](#nspin) | [nspinor](#nspinor) | [syms](#syms)| [sys](#sys) 

### C. Properties details

#### rhor

```
Type       : Double
Size       : (n1, n2, n3)
Description: Charge density on the real-space FFT grid
```

#### Vxc

```
Type       : Double
Size       : (nb, nk, nspin)
Description: Exchange-correlation matrix elements / values used in QP correction
Note       : Units follow the loader path; qp.fout treats Vxc as eV in the table
```

#### ev

```
Type       : Double
Size       : (nb, nk, nspin)  (shape may vary slightly by loader)
Description: Kohn–Sham eigenvalues
Note       : QE path converts to Rydberg internally where applicable
```

#### psig

```
Type       : Cell
Size       : {nk, nspin}, each entry complex (ng, nb)
Description: Wavefunctions in reciprocal space (G basis per k)
```

#### occupation

```
Type       : Double
Size       : (nb, nk, nspin)
Description: Band occupation numbers
```

#### reciprocal_grid_info

```
Type       : Struct
Description: FFT / G-mapping for wavefunctions
Typical fields:
  fftgrid, vol, idxnz (cell per k), wfncut, xyz
```

#### nkibz / kibz / kweight

```
Type       : Integer / Double / Double
Description: Number of IBZ k-points, coordinates, and weights
```

#### nspin / nspinor

```
Type       : Integer
Description: Spin channels and spinor components
```

#### syms

```
Type       : Struct
Description: Symmetry information imported from ground state
Typical fields:
  nsym, nrot, is_t_rev, mtrx, indsub, kgzero, ...
```

#### sys

```
Type       : Struct
Description: Cell / grid metadata, only used in ISDF with `exxmethod'=`kmeans'
Typical fields:
  ng, nr, ne, n1, n2, n3, supercell, qk, vol
```


### D. Usage

```matlab
data = load(.../SAVE/data.mat', 'data').data;
% After input_driver / service_driver, prefer service getters for kernels.
```

---

---

## Service layer (`service/+*/+base/*_m`)

Service objects are MATLAB `classdef` types held by persistent managers.  
They are constructed in `service_driver` and (partially) cached in `relay_stage.mat`.

Typical access patterns:

```matlab
sys  = system.get();
wf   = wave_functions.get();
coul = coulomb.get();
rlat = lattice.manager('r_lat');
dlat = lattice.manager('d_lat');
kobj = lattice.manager('k');
qobj = lattice.manager('q');
```

### Cache policy (relay)

**Cached in `relay_stage.mat`:**  
`system`, `symmetry`, `FFT`, `wave_functions` (and related: `pair_symmetry`, `timing` when enabled)

**Rebuilt every run from current `config`:**  
`lattice`, `coulomb`, ISDF

---

### `system.base.system_m`

#### A. Introduction

Electronic / QP bookkeeping: eigenvalues, Vxc, occupations, and degeneracy tables.

#### B. Properties list

[nb](#system-nb) | [nk](#system-nk) | [nspin](#system-nspin) | [Eo](#eo) | [Vxc](#system-vxc) | [first_index_in_degeneracy](#first_index_in_degeneracy) | [num_index_in_degeneracy](#num_index_in_degeneracy) | [degeneracy_indices_len](#degeneracy_indices_len)

#### C. Properties details

#### nb / nk / nspin {#system-nb}

```
Type       : Integer scalars
Description: Number of bands, IBZ k-points, spin channels
```

#### Eo

```
Type       : Double
Size       : (nb, nk, nspin)
Description: Mean-field / KS energies stored by the system service
```

#### Vxc {#system-vxc}

```
Type       : Double
Size       : (nb, nk, nspin)
Description: Exchange-correlation contribution used in Eqp0 = Eo + Sig - Vxc
```

<!-- #### Eqp {#system-eqp} -->
<!--  -->
<!-- ``` -->
<!-- Type       : Double -->
<!-- Size       : (nb, nk, nspin) -->
<!-- Description: Working buffer for quasiparticle energies inside system service -->
<!-- ``` -->

#### f

```
Type       : Double
Size       : (nb, nk, nspin)
Description: Occupations
```

<!-- #### qptype -->
<!--  -->
<!-- ``` -->
<!-- Type       : String -->
<!-- Description: QP type tag (default "HF") -->
<!-- ``` -->

#### first_index_in_degeneracy / num_index_in_degeneracy / degeneracy_indices_len

```
Type       : Cell / Cell / int32
Description: Degeneracy grouping per k-point for averaging / output
```

---

### `wave_functions.base.wf_m`

#### A. Introduction

Real-space wavefunction container.

#### B. Properties

Dimension counts (`nb` / `nk` / `nspin` / `n_spinor` / `ng`) are copies owned by
other modules; stored on `wf_m` so kernels can size `c` without cross-service lookups.

| Name | Type | Size | Description |
|------|------|------|-------------|
| `c` | complex double | `(nc, nb, nk, nspin)` | Wavefunction coefficients on the FFT grid (`space = "r"` by default) |
| `nc` | int | scalar | Real-space FFT points `n1*n2*n3 = prod(FFT.fftgrid)` |
| `space` | string | scalar | Representation flag; default `"r"` |
| `nb` | int | scalar | Number of bands → `system.nb` |
| `nk` | int | scalar | IBZ k-points → `lattice k.nibz` |
| `nspin` | int | scalar | Spin channels → `system.nspin` |
| `n_spinor` | int | scalar | Spinor components (1 / 2) → `data.nspinor` (ctor sets 1; spinor path not fully wired) |
| `ng` | int | scalar | G-vector count → `r_lattice.ng` (unused while `space = "r"`; stays 0 after driver) |

---

### `lattice.base.d_lattice_m` / `r_lattice_m` / `bz_samp_m`

#### A. Introduction

Lattice and Brillouin-zone sampling. Accessed via `lattice.manager('d_lat'|'r_lat'|'k'|'q')`.

#### B. `d_lattice_m` (direct lattice)

| Name | Type | Size | Description |
|------|------|------|-------------|
| `a1a2a3` | double | `(3, 3)` | Lattice vectors as columns |
| `DL_vol` | double | scalar | Cell volume |
| `atom_pos` | double | `(:,:,3)` | Atomic positions |
| `atom_symbol` | cell | — | Element symbols |

#### C. `r_lattice_m` (reciprocal lattice)

| Name | Type | Size | Description |
|------|------|------|-------------|
| `b1b2b3` | double | `(3, 3)` | Reciprocal vectors as columns |
| `RL_vol` | double | scalar | Reciprocal-cell volume |
| `ng` | int | scalar | Number of G vectors in cutoff |
| `Ggrid_RLU` | int | `(ng, 3)` | G in reciprocal lattice units |
| `Ggrid_Cart` | double | `(ng, 3)` | G in Cartesian |
| `idxnz` | int | `(ng, 1)` | Map to full FFT grid |
| `G_rot` | int | `(ng, nsym)` | Symmetry action on G |
| `qindx_S` / `qindx_X` / `qindx_C` / `qindx_B` | int | `(nibz/nbz, nbz, 2)` | Momentum-transfer index tables |
| `n_g_shell` | int | scalar | Number of G shells |
| `num_index_in_each_shell` | int | `(n_g_shell,)` | Shell sizes |
| `first_index_in_each_shell` | int | `(n_g_shell,)` | Shell start indices |

#### D. `bz_samp_m` (k / q sampling)

| Name | Type | Size | Description |
|------|------|------|-------------|
| `nibz` | int | scalar | Number of IBZ k-points |
| `nbz` | int | scalar | Number of full-BZ k-points |
| `kpt_RLU` | double | `(nibz, 3)` | IBZ coordinates (RLU) |
| `kpt_Cart` | double | `(nibz, 3)` | IBZ coordinates (Cartesian) |
| `kptbz_RLU` | double | `(nbz, 3)` | Full-BZ coordinates (RLU) |
| `kptbz_Cart` | double | `(nbz, 3)` | Full-BZ coordinates (Cartesian) |
| `weights` | double | `(nibz,)` | IBZ integration weights |
| `nstar` | int | `(nibz,)` | Star size per IBZ point |
| `star` | int | `(nibz, nsym)` | Symmetry indices in each star |
| `sstar` | int | — | Star membership map |
| `s_table` | int | `(nibz, nsym)` | `(ik_ibz, isym) → ik_bz` |
| `k_table` | int | `(nibz, nsym)` | `(ik_ibz, isym) → ik_bz` (alt) |
| `bz2ibz` | int | `(nbz,)` | BZ → IBZ index |
| `bz2rot` | int | `(nbz,)` | BZ → rotation index |
| `ibz2bz` | int | — | IBZ → representative BZ index |

---

### `coulomb.base.coulomb_m`

#### A. Introduction

Bare / truncated Coulomb kernel in reciprocal space.

#### B. Properties

| Name | Type | Size | Description |
|------|------|------|-------------|
| `trunc_method` | int | scalar | Truncation scheme (see CUTOFFS) |
| `trunc_param` | double | scalar | Truncation parameter |
| `coulomb_ng` | int | scalar | Number of G vectors used for Coulomb |
| `bare_qpg` | double | `(ng, nqibz)` | \|q+G\| |
| `vcoul` | double | `(ng, nqibz)` | Coulomb kernel values |
| `vcoul0` | double | scalar | Coulomb at q = G = 0 |

---

### `FFT.base.FFT_m`

#### A. Introduction

Real-space FFT grid and symmetry permutation tables.

#### B. Properties

| Name | Type | Size | Description |
|------|------|------|-------------|
| `fftgrid` | int | `(1, 3)` | Real-space FFT dimensions `(n1, n2, n3)` |
| `nr` | int | scalar | `n1*n2*n3` |
| `Rgrid_RLU` | double | `(nr, 3)` | R-point coordinates in lattice units |
| `R_rot` | int | `(nr, nsym)` | Symmetry permutation of R points |
| `R_rot_inv` | int | `(nr, nsym)` | Inverse permutation |
| `nGo` | int | scalar | G0-shell count for FFT box |
| `G_table` | — | — | G0-shell helpers for FFT box |

---

### `symmetry.base.symm_m`

#### A. Introduction

Crystal symmetries (rotations / time-reversal flags).

#### B. Properties

| Name | Type | Size | Description |
|------|------|------|-------------|
| `nsym` | int | scalar | Number of symmetries used |
| `nrot` | int | scalar | Number of rotations |
| `is_t_rev` | logical  | scalar  | Time-reversal flag |
| `rot_mtrx_RLU_G` | double | `(3, 3, nsym)` | Rotation in RLU (G convention) |
| `rot_mtrx_RLU_R` | double | `(3, 3, nsym)` | Rotation in RLU (R convention) |
| `rot_mtrx_Cart` | double | `(3, 3, nsym)` | Rotation in Cartesian |
| `inv_rot_index` | int | `(nsym, 1)` | `inv_rot_index(i1) = i2` such that `S(i1) * S(i2) = I` (i.e. `S(i2) = S(i1)^{-1}`) |

---

### `isdf.base.isdf_m`

#### A. Introduction

ISDF pool entry (one object per channel / id, e.g. `vc`, `vn`, `nn`).  
Access via the ISDF package API after `isdf.driver(data, config)`.

#### B. Properties

| Name | Type | Size | Description |
|------|------|------|-------------|
| `id` | int | scalar | Pool identifier |
| `desc` | string | scalar | Channel tag (`vc` / `vn` / `nn`, …) |
| `nisdf` | int | scalar | Number of interpolation points |
| `nrange1` | int | `(2,)` or range | Band range for first product index |
| `nrange2` | int | `(2,)` or range | Band range for second product index |
| `coeff_seper` | complex / double | `(nisdf, nb, nk, nspin)` | Separation coefficients |
| `tildeVq` | complex / double | `(nisdf, nisdf, …)` | Compressed Coulomb / interaction |
| `helperqG` | — | — | Helper functions in G space |
| `CCHq` | — | — | Overlap metadata |
| `CCHq_trunc_factors` | — | — | Truncation metadata for `CCHq` |
| `fftgrid_c` | int | `(1, 3)` | Coarse FFT grid |
| `R_rot_coarse` | int | — | Coarse-grid R rotation |
| `R_rot_extra` | int | — | Extra / adaptive R rotation |
| `N_coarse` | int | scalar | Adaptive coarse count |
| `N_extra` | int | scalar | Adaptive extra count (`N_coarse + N_extra = nisdf`) |

---


---

## `E` (return value of `qp.launcher`)

### A. Introduction

`E` is the lightweight result struct returned by `qp.launcher` / `qp_driver`.  
Human-readable tables are written by `qp.fout` to `qp.dat`.

**Internal energies in `E.Ex` / `E.Esx_x` / `E.Ecoh` / `E.Eqp` are in Rydberg.**  
`qp.fout` converts them to **eV** for the text table (`ry2ev = 13.60569253`).

### B. Properties

| Name | Type | Size | Description |
|------|------|------|-------------|
| `Eqp` | double | `(nband_window, nk)` | `Ex + Esx_x + Ecoh` for selected bands (Ry) |
| `Ex` | double | `(nband_window, nk)` | Exact exchange Σ_x (Ry) |
| `Esx_x` | double | `(nband_window, nk)` | Screened-exchange correction beyond X (COHSEX SX-X, or full-frequency residual) |
| `Ecoh` | double | `(nband_window, nk)` | Coulomb-hole / integral contribution (Ry) |
| `Eqp0` | double / `[]` | — | Not filled by `qp.launcher`; `qp.fout` computes `Eqp0 = Eo + (Ex+Esx_x+Ecoh) - Vxc` for the table |
| `fout` | string / char | scalar | Path to written `qp.dat` (empty if write failed) |

### C. Related helpers

```matlab
E    = qp.launcher(config);
fout = qp.fout(E, config);       % also called inside launcher
qp.summary(0, config);           % pre
qp.summary(1, config, E, t);     % post
```

---

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-08 | Zhengbang | Initialize for data + service `*_m` + `E` |
