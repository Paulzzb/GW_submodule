# benchmarks/

No QE wavefunctions are shipped. Run GS yourself, then GW.

| Dir | What | ecut | GW |
|-----|------|-----:|----|
| `Si8/` | Si 8-atom | 40 | full-freq + ISDF/direct |
| `Si64/` | Si 2×2×2 | 40 | COHSEX + ISDF/direct |
| `STO3/` | SrTiO₃ 5-atom | 80 | full-freq + ISDF/direct |
| `STO3_8/` | SrTiO₃ 2×2×2 | 80 | COHSEX + ISDF/direct |

UPFs: Si in `Si64/`, STO in `STO3_8/`. Each case: `scf.in` → `nscf.in` → `pp_in`, plus `test_isdf` / `test_dir`.

## Run (SLURM)

```bash
mkdir -p log
sbatch slurm_qe_gs_small.sh    # Si8, STO3
sbatch slurm_qe_gs_large.sh    # Si64, STO3_8
sbatch slurm_gw_small.sh       # default: test_isdf → ./isdf/qp.dat
sbatch slurm_gw_large.sh
NAMELIST=test_dir sbatch slurm_gw_small.sh   # dense → ./direct/qp.dat
```

Overrides: `CASES=…`, `NPROC=…`, `PW=…`, `MATLAB_BIN=…`.

## Run (manual)

```bash
cd benchmarks/Si8
pw.x < scf.in > scf.out && pw.x < nscf.in > nscf.out && pw2bgw.x < pp_in > pp.out
```

```matlab
QPstartup
cd benchmarks/Si8
input_driver('./test_isdf'); load('./SAVE/config.mat','config'); E = qp.launcher(config);
```

`STO3` and `STO3_8` both use QE prefix `STO3` (each keeps its own `STO3.save`). Large / dense jobs need lots of memory.
