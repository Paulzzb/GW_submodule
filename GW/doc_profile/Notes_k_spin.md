# TEST

We generate some new class, and add some new properties to old classes to support the calculation.

In the calculation, there are reducible and irreducible BZs, the latter one is denoted by IBZ.
We need to sample the IBZ to get the finite set of k-points, named as reducible
k-points set.
However, due to the symmetry, this set can be further reduced to a minimal set, which is named as irreducible k-points set.

Normally, the relationship is
\[
kbz(i) + G(iGolist(i)) = \hat{S}_(ind_rotation(i)) * kibz(ind_kbz(i)), 
\]
All of these are stored in the class `@bz_sampling`.

When there goes to the momentum index $\qq$, we need to calculate $\kk' = \kk - \qq$, and the relation is as followed:
\[
k_ibz(i) - kbz(j) = kbz(qindx_S(i, j, 1)) + G(qindx_S(i, j, 2))
\]
The mapping qindx_S is also stored in the class `@bz_sampling`.

When we calculate the densities, we are indeed facing
...

## @bz_sampling

### A. Class introduction

`@bz_sampling` is a supporting class designed to manage Brillouin zone (BZ)
related information.
Irreducible Brillouin zone (IBZ) sampling for spin- and momentum-resolved quasiparticle calculations and all related mapping of k-points' indices are
included in this class.  
It is constructed during the `input_driver` stage and passed to downstream computational routines to manage $\mathbf{k}$- and $\mathbf{q}$-space interpolation and symmetries.  
It enables efficient lookup of index mappings and interpolation tables necessary for evaluating quantities like the self-energy $\Sigma(\mathbf{k}, \omega)$.

### B. Properties list

[nibz](#nibz) |
[nbz](#nbz) |
[kpt](#kpt) |
[kptbz](#kptbz) |
[kptweights](#kptweights) |
[bmatrix](#bmatrix) |
[kbz2kibz_ind_rotation](#kbz2kibz_ind_rotation) |
[kbz2kibz_ind_kbz](#kbz2kibz_ind_kbz) |
[kbz2kibz_ind_Go](#kbz2kibz_ind_go) |
[nstar](#nstar) |
[star](#star) |
[sstar](#sstar) |
[s_table](#s_table) |
[k_table](#k_table) |
[qindx_S](#qindx_s) |
[qindx_C](#qindx_c) |
[qindx_X](#qindx_x) |
[qindx_B](#qindx_b) |
[nGo](#ngo) |
[iGolist](#igolist)


---

### C. Properties details


#### nibz

```
Type       : integer
Size       : scaler
Description: Number of irreducible k-points
```

#### nbz

```
Type       : integer
Size       : scaler
Description: Number of reducible k-points
```

#### kpt

```
Type       : double precision
Size       : (nibz, 3)
Description: irreducible k-points set, in fractional coordinates.
```

#### kptbz

```
Type       : double precision
Size       : (nbz, 3)
Description: reducible k-points set, in fractional coordinates.
```

#### kptweights

```
Type       : double precision
Size       : (nibz, 1)
Description: weights for each irreducible k-points.
```

#### bmatrix

```
Type       : double precision
Size       : (3, 3)
Description: reciprocal lattice vectors.
```

#### kbz2kibz_ind*

```
Type       : double precision
Size       : (3, 3)
Description: reciprocal lattice vectors.
```


## GWinfo

