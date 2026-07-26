# +isdftest — SVD-based ISDF test package

Copy of `+isdf` with a separate manager pool and modified `tildeVq` construction.

## Difference from +isdf

| Step | +isdf | +isdftest |
|------|-------|-----------|
| Helper in R | `MCHq / CCHq` | `MCHq * C^{-1/2}` via `prod_C_inv_sqrt('prod','r',...)` |
| `tildeVq` | `helperqG' V helperqG` | same form, different helper |
| Application | `c' * tildeVq * c` | `get_rho_xalpha` returns $c_t=C^{-1/2}c$; then `c_t' tildeVq c_t` |
| Per-q cache | — | `CCHq`, `CCHq_inv_sqrt(:,:,iq)` on `isdftest_m` |

## Key entry points

- `isdftest.gen_tildeVq(id)` — build; optional `s_cut` truncates small eigenvalues of `CCHq`
- `isdftest.get_rho_xalpha(id, param)` — returns $c_t=C^{-1/2}c$ (requires prior `gen_tildeVq`)
- `isdftest.contract_tildeVq(id, iqibz, c)` — optional wrapper for `c' * tildeVq * c`
- `isdftest.prod_C_inv_sqrt('set', CCHq, s_cut)` / `'prod', 'l'|'r', A` — apply truncated $C^{-1/2}$ without forming the full matrix

See `GW/test_profile/demo_isdf_algorithm.tex` and `demo.m` for algorithm notes and usage.
