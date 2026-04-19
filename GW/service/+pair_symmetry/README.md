# pair_symmetry Module

This module stores and constructs symmetry mapping data for composite k-point pairs $(k_1, k_2)$ on the full BZ mesh.

The current implementation is centered on the class `pair_symmetry.base.pair_symm_m` and the builder `pair_symmetry.driver`.

## 1. Purpose

For each target pair of BZ indices

$$
(i_{k_1'}, i_{k_2'}) \, ,
$$

the module tries to identify a representative/source pair

$$
(i_{k_1}, i_{k_2})
$$

and a symmetry operation $S$ such that

$$
k_1' = S k_1 + G_0,
$$

$$
k_2' = S k_2 + G_0,
$$

with the **same** reciprocal-lattice shift $G_0$ for both components.

This shared-$G_0$ condition is the core structural constraint of the module.

## 2. Main Data Object

The module stores one object of class `pair_symmetry.base.pair_symm_m` in the manager cache.

Main properties:

- `assigned`: whether the object has been filled with meaningful mapping data.
- `allocated`: whether storage has been allocated.
- `nk`: number of full-BZ k-points.
- `representation(irep, :)`: representative/source pair table.
- `weights(irep)`: normalized weight of each representative.
- `mapping(ik1p, ik2p, :)`: dense lookup table for every target pair.
- `g0_table(iG0, :)`: table of unique shared reciprocal shifts.
- `flagconj(ik1p, ik2p)`: whether the target pair was filled through the swapped-order branch.

## 3. representation

`representation` has size

$$
N_{\mathrm{rep}} \times 2.
$$

Each row stores one representative/source pair:

$$
\texttt{representation}(i_{\mathrm{rep}}, :) = [i_{k_1}, i_{k_2}].
$$

These are the orbit seeds chosen by the driver while scanning the full pair grid.

## 4. weights

`weights` has size

$$
N_{\mathrm{rep}} \times 1.
$$

For representative index $i_{\mathrm{rep}}$, define orbit coverage count

$$
c_{i_{\mathrm{rep}}} = \#\{(i_{k_1'}, i_{k_2'}) : \texttt{mapping}(i_{k_1'}, i_{k_2'}, 1) = i_{\mathrm{rep}}\}.
$$

The module stores the normalized weight

$$
w_{i_{\mathrm{rep}}} = \frac{c_{i_{\mathrm{rep}}}}{n_k^2}.
$$

So `weights` represents how much of the full pair grid is represented by each row in `representation`.

## 5. mapping

`mapping` has size

$$
n_k \times n_k \times 3.
$$

For each target pair `(ik1p, ik2p)`:

- `mapping(ik1p, ik2p, 1)` = `irep`
- `mapping(ik1p, ik2p, 2)` = `isym`
- `mapping(ik1p, ik2p, 3)` = `iG0`

Interpretation:

1. Read the representative pair from `representation(irep, :)`.
2. Read the symmetry operation index `isym`.
3. Read the reciprocal shift from `g0_table(iG0, :)`.

Then reconstruct the target pair by

$$
k_1' = S_{i_{\mathrm{sym}}} k_{1,\mathrm{rep}} + G_0,
$$

$$
k_2' = S_{i_{\mathrm{sym}}} k_{2,\mathrm{rep}} + G_0.
$$

So the first channel of `mapping` is **not** a k-point index. It is an index into `representation`.

## 6. g0_table

`g0_table` stores unique integer reciprocal-lattice shifts in RLU form:

$$
G_0 = (g_1, g_2, g_3), \qquad g_i \in \mathbb{Z}.
$$

If

$$
i_{G_0} = \texttt{mapping}(i_{k_1'}, i_{k_2'}, 3),
$$

then

$$
G_0 = \texttt{g0\_table}(i_{G_0}, :).
$$

This indirection keeps `mapping` compact and purely integer-valued.

## 7. flagconj

`flagconj(ik1p, ik2p)` is logical.

It is set to `true` when the target pair is filled through the swapped-order branch

$$
(k_1', k_2') \leftrightarrow (k_2', k_1'),
$$

rather than the direct branch.

Current practical meaning:

- `false`: direct orbit fill
- `true`: swapped target fill

This is intended to track the branch associated with pair-order exchange, which is typically conjugate-related in later formulas.

## 8. Driver Construction Logic

The main builder is `pair_symmetry.driver`.

Its internal algorithm is:

1. Read full-BZ k-points from `lattice.manager('k', 'get')`.
2. Read symmetry matrices from `symmetry.manager('get')`.
3. Scan all target pairs `(ik1, ik2)`.
4. If the pair is still unassigned, choose it as a new representative.
5. Expand its orbit over all symmetry operations.
6. For each reachable pair, store:
   - representative index `irep`
   - symmetry index `isym`
   - shared shift index `iG0`
7. Also fill the swapped target and mark `flagconj = true` there when appropriate.
8. Count `mapping(:, :, 1)` and normalize by `nk^2` to generate `weights`.

## 9. Shared-G0 Matching Formula

For one representative pair $(k_1, k_2)$ and one symmetry $S$, the code computes

$$
S k_1, \qquad S k_2.
$$

It then searches BZ mesh indices `(ia, ib)` such that

$$
k(ia) - S k_1 = G_0,
$$

$$
k(ib) - S k_2 = G_0,
$$

with the same integer vector $G_0$.

Numerically, this is checked through

$$
G_{0,1} = \mathrm{round}(k(ia) - S k_1),
$$

$$
G_{0,2} = \mathrm{round}(k(ib) - S k_2),
$$

and acceptance requires

$$
\| (k(ia) - S k_1) - G_{0,1} \| \le \mathrm{tol},
$$

$$
\| (k(ib) - S k_2) - G_{0,2} \| \le \mathrm{tol},
$$

$$
\| G_{0,1} - G_{0,2} \| \le \mathrm{tol}.
$$

Among valid matches, the implementation prefers:

1. $G_0 = 0$
2. otherwise the smallest $\|G_0\|_1$

## 10. Current Files

Current module files:

- `+base/pair_symm_m.m`: storage class
- `driver.m`: builds and stores mapping data
- `manager.m`: persistent cache manager
- `get.m`: wrapper for `manager('get')`
- `save2mod.m`: wrapper for `manager('save2mod', ...)`
- `free.m`: wrapper for `manager('free')`
- `save_k1k2_result.m`: MAT-file saving helper
- `validate_k1k2_mapping.m`: transitional validator

## 11. Transitional Note

`validate_k1k2_mapping.m` is currently a transitional checker and does not yet fully consume the newer `mapping(:, :, 1:3) + g0_table + flagconj` contract in its most natural form.

So at this stage:

- the storage model is the new one,
- the builder follows the new one,
- the validator is only partially updated.

This should be kept in mind when extending the module further.
