# Si_gamma verification groups

Shared groundstate: `Si_gamma/qe.save` (siblings use `groundstate_dir = '../Si_gamma/qe.save'`).

| Group | Purpose | `frequency_dependence` | Cases |
|------|---------|--------------------------|--------|
| G1 | Three ISDF methods (`exxmethod`) | `2` | `Si_gamma` (pseudo), `Si_gamma_qrcp`, `Si_gamma_kmeans` |
| G2 | ISDF vs dense | `2` | `Si_gamma` (`isisdf=1`), `Si_gamma_dir` (`isisdf=0`) |
| G3a | Cauchy on/off (full-frequency) | `2` | `Si_gamma_ff_cauchy`, `Si_gamma_ff_nocauchy` |
| G3b | Cauchy on/off (COHSEX Gamma) | `-2` | `Si_gamma_cohsex_cauchy`, `Si_gamma_cohsex_nocauchy` |
| G4 | Ratio × exact_ch / FF (`exxmethod=pseudo`) | see below | `Si_gamma_16*` |

## G4 — ratio / exact_ch / FF

| Case | ratios (vc/vn/nn) | `exact_ch` | `fdep` | Checkpoint reuse |
|------|-------------------|------------|--------|------------------|
| `Si_gamma_162416_exact` | 16 / **24** / 16 | `.true.` | `-2` | vc,nn ← hub; **builds vn** |
| `Si_gamma_162416_sum` | 16 / 24 / 16 | `.false.` | `-2` | vc,nn ← hub; vn ← `162416_exact` |
| `Si_gamma_162416_ff` | 16 / 24 / 16 | (n/a FF) | `2` | same as sum |
| `Si_gamma_161616_exact` | 16 / 16 / 16 | `.true.` | `-2` | all adaptive ckpts ← hub |
| `Si_gamma_161616_sum` | 16 / 16 / 16 | `.false.` | `-2` | all ← hub |
| `Si_gamma_161616_ff` | 16 / 16 / 16 | (n/a FF) | `2` | all ← hub |

Hub `Si_gamma` is 16/16/16 (+ FF). Run order: hub → … → **`162416_exact` before other 162416\*** → 161616\*.

Links: `example/link_isdf_checkpoints.sh` (also called from `run_all` after each clean).

Run:

```matlab
cd <repo>/example
run_all
collect_qp_energies
```

Or SLURM:

```bash
mkdir -p log
sbatch run_all_shared_save.sh
```
