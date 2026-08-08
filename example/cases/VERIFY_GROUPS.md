# Si_gamma verification groups

Shared groundstate: `Si_gamma/qe.save` (siblings use `groundstate_dir = '../Si_gamma/qe.save'`).

| Group | Purpose | `frequency_dependence` | Cases |
|------|---------|--------------------------|--------|
| G1 | Three ISDF methods (`exxmethod`) | `2` | `Si_gamma` (pseudo), `Si_gamma_qrcp`, `Si_gamma_kmeans` |
| G2 | ISDF vs dense | `2` | `Si_gamma` (`isisdf=1`), `Si_gamma_dir` (`isisdf=0`) |
| G3a | Cauchy on/off (full-frequency) | `2` | `Si_gamma_ff_cauchy`, `Si_gamma_ff_nocauchy` |
| G3b | Cauchy on/off (COHSEX Gamma) | `-2` | `Si_gamma_cohsex_cauchy`, `Si_gamma_cohsex_nocauchy` |

Run all listed cases:

```matlab
cd <repo>/example
run_all
collect_qp_energies
```
