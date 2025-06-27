# Introduction and reference of classes in QP source code 

We introduce three class in source code, named @GWinfo, @gvec, and @QPenergy.
They are used to

- @GWinfo: save groundstate informations.
- @gvec: save reciprocal space information.
- @QPenergy: save quasiparticle energies' information

For the convenience of developers to query, we summarize the important information as follows

---
---

## @GWinfo

### A. Class introduction

`@GWinfo` is a core data structure in the software package that stores all ground-state information required for quasiparticle (QP) calculations.
Its primary function is to serve as an information carrier, transferring ground-state data between `input_driver` and `qp_driver`,
and also acting as a key interface within `qp_driver`.

Parameter Convention

`Ncoul`: Number of reciprocal space grid points under the given Coulomb cutoff.

### B. Properties list

[coulG](#coulg) |
[coulG0](#coulg0) |
[supercell](#supercell) |
[bdot](#bdot) |
[qk](#qk) |
[vol](#vol) |
[ne](#ne) |
[gvec](#gvec) |
[gvec2](#gvec2) |
[gvecrho](#gvecrho) |
[Vxc](#vxc) |
[rho](#rho) |
[idxnz](#idxnz) |
[ev](#ev) |
[psig](#psig) |
[psir](#psir) |
[Ggrid4psig](#ggrid4psig) |
[occupation](#occupation) |

---
### C. Properties details

#### coulG

```
Type       : Double
Size       : (Ncoul, 1)
Description: Coulomb kernel under specific approximation method
and specific truncation in reciprocal space
```

#### coulG0
```
Type       : Double 
Size       : scalar
Description: Coulomb potential at q=G=0
```

#### supercell
```
Type       : Double
Size       : (3, 3)
Description: Supercell lattice matrix in real space, C = [a1, a2, a3]'
```

#### bdot (Deprivated, not initialized)
```
Type       : Double
Size       : (3, 3)
Description: A matrix to perform vector dot in reciprocial space.
Note       : Suppose B = [b1, b2, b3] is reciprocal space
primitive vectors, then bdot = B'*B.
```

#### qk (Deprivated, not initialized)
```
> Todo
Type       : Double 
Size       : ...
Description: 
Note       : 
```

#### vol
```
Type       : Double
Size       : scalar
Description: Volume of the supercell
```

#### ne (Deprivated, not initialized)
```
Type       : integer
Size       : scalar
Description: Total number of electrons
```

#### gvec
```
Type       : @gvec 
Size       : Not applicable 
Description: Reciprocal lattice vector information, see (#@gvec) for detail
Note       : Used in expansion of wavefunctions
```

#### gvec2 (Deprivated, not initialized)
```
Type       : @gvec
Size       : 
Description: Was used on description of wavefunctions grid in reciprocal space.
Note       : 
```

#### gvecrho
```
Type       : @gvec
Size       : Unknown
Description: reciprocal vectors grid, used for charge density (typically finer grid)
Note       : Different cutoff from wavefunction basis
```

#### Vxc
```
Type       : Double 
Size       : (Number of required bands, 1)
Description: Exchange-correlation potential in real space
Note       : Imported from ground-state solver output
```

#### rho
```
Type       : Double Complex
Size       : (gvecrho.ng, 1)
Description: Charge density in reciprocal space, always used with gvecrho
Note       : 
```

#### idxnz (Deprivated, not initialized)
```
Type       : Integer 
Size       : 
Description: 
Note       : 
```

#### ev
```
Type       : Double 
Size       : (nband, 1)
Description: ground-state band energies
```

#### Ggrid4psig
```
Type       : struct 
Size       : Not appliable  
Description: Grid mapping information between truncated G-vectors and real-space grid
Note       : Use get_wavefunc_real(psig, Ggrid4psig) to get psir
```

#### psig
```
Type       : Double Complex
Size       : Not important (hhh) 
Description: Wavefunctions in reciprocal space 
Note       : Use get_wavefunc_real(psig, Ggrid4psig) to get psir
```

#### psir

```
Type       : complex array
Size       :  
Description: Wavefunctions in real space
Note       : Use get_wavefunc_real(psig, Ggrid4psig) to get psir.
Empty when outputed from input_driver, initialized in qp_driver.
```

#### occupation
```
Type       : Double 
Size       : (nband, 1)
Description: Band occupation numbers (0 to 1)
Note       : Determines which bands are treated as occupied/unoccupied
```

---

### C. Functions and routines

Currently, @GWinfo is a 'storage' room.
Please just use GWinfo.property_name to extract information.

---
---

## @gvec

### A. Class Introduction

`@gvec` is a class that contains information about the truncated reciprocal-space grids.
Developers typically invoke it to perform a forward or backward FFT through
a set of functions provided in the folder `QP_root/common`.

### B. Properties List

[ng](#ng) |  
[fftgrid](#fftgrid) |  
[nfftgridpts](#nfftgridpts) |  
[idxnz](#idxnz) |  
[components](#components) |

---

### C. Properties Details

#### ng

```
Type       : Integer
Size       : Scalar
Description: Number of reciprocal grid vectors remaining after applying the energy cutoff
```

#### fftgrid

```
Type       : Integer
Size       : (3, 1)
Description: Inital real-space grid size in each direction (n1, n2, n3)
```

#### nfftgridpts
```
Type       : Integer
Size       : Scalar
Description: Total number of FFT grid points, e.g., n1*n2*n3
```

#### idxnz
```
Type       : Integer Array
Size       : (ng, 1)
Description: Linear index mapping of retained grid points (1:ng) to the full FFT grid (n1, n2, n3).
Note: We use MATLAB convention
```
> Todo: implementing detail

#### components
```
Type       : Integer
Size       : (ng, 3)
Description: Component-wise mapping of grid points in 3D space, i.e., (i1, i2, i3) with i1 ∈ [1,n1] etc.
```

---

### D. Functions and Routines

To construct a `@gvec` object:

```matlab
gvec_obj = gvec(gvecinput)
```

Where `gvecinput` must contain:

- `n1, n2, n3`: FFT grid dimensions
- `ecut`: Energy cutoff for grid point truncation
- `supercell`: Lattice matrix in real space

---
---

## @QPenergy

### A. Class Introduction

`@QPenergy` is the primary data structure that stores and organizes all quasiparticle (QP) energy-related outputs during and after the quasiparticle calculation.  
It is initialized in `qp_driver`, and acts as a central container passed through `qp_driver` and post-processing stage.
Post-processing functions like `getEqp`, `shiftenergy`, and `GWfout` and contained in this class.
**All the energies have unit ev**

---

### B. Properties List

[bandindices](#bandindices) |
[freq_dep](#freq_dep) |
[freq_dep_method](#freq_dep_method) |
[ev](#ev) |
[Ex](#ex) |
[Esex_x](#esex_x) |
[Ecoh](#ecoh) |
[Eres](#eres) |
[Eint](#eint) |
[Sig](#sig) |
[Vxc](#vxc) |
[Eqp](#eqp) |
[TOL_DEGENERACY](#tol_degeneracy) |
[fout](#fout)

---

### C. Properties Details

#### bandindices

```
Type       : Integer Array
Size       : (nbands_out)
Description: Indices of bands involved in QP calculation,
parsed by parameters energy_band_index_min and energy_band_index_max from user-provide input
```                   

#### freq_dep

```
Type       : Integer
Size       : Scalar
Description: Frequency-dependence type (0=static, 1=GPP, 2=full-frequency)
```

#### freq_dep_method

```
Type       : Integer
Size       : Scalar
Description: Method used for frequency grid generation
```

#### ev

```
Type       : Double
Size       : (nband_out, 1)
Description: DFT eigenvalues of selected bands
```

#### Ex

```
Type       : Double
Size       : (nband_out, 1)
Description: Exchange energy for each band
```

#### Esex_x

```
Type       : Double
Size       : (nband_out, 1)
Description: Static screened exchange component (freq_dep = 0 or 1)
```

#### Ecoh

```
Type       : Double
Size       : (nband_out, 1)
Description: Coulomb hole contribution (freq_dep = 0 or 1)
```

#### Eres

```
Type       : Double
Size       : (nband_out, 1)
Description: Residue term from full-frequency GW (freq_dep = 2)
```

#### Eint

```
Type       : Double
Size       : (nband_out, 1)
Description: Integral term from full-frequency GW (freq_dep = 2)
```

#### Sig

```
Type       : Struct or Cell
Size       : varies
Description: Self-energy function evaluated at selected bands
```

#### Vxc

```
Type       : Double
Size       : (nband_out, 1)
Description: Exchange-correlation energies from ground-state calculation
```

#### Eqp

```
Type       : Double
Size       : (nband_out, 1)
Description: Final quasiparticle energy bands (after correction)
```

#### TOL_DEGENERACY

```
Type       : Double
Size       : Scalar (constant)
Description: Tolerance for detecting energy degeneracies, on default 1e-6*13.60569
```

#### fout

```
Type       : String
Size       : Scalar
Description: Default file name to write QP energy output
```

---
---



### 4. Changing log


| Date | Name   | Changes       |
|------------|--------|------------|
| 2025-06-27 | Zhengbang | Initialize|
---
