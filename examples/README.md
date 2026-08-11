# examples/ — public ISDF vs dense demo

Small Gamma-Si showcase: same groundstate, compare **ISDF** vs **dense** for
full-frequency (`fdep=2`) and COHSEX (`fdep=-2`). Timing and accuracy are
the point — each case keeps its own `SAVE/` so setup cost is not shared away.

User guide: [`../doc_user_EN/README.md`](../doc_user_EN/README.md) · [`../doc_user_ZH/README.md`](../doc_user_ZH/README.md).  
Internal verification suites live under `tests/` (not this tree).

## Layout

```text
examples/
  run_all.m / run_all_profiled.m / run_batch.sh
  run_cohsex.m
  cases/
    qe.save/              % shared QE groundstate only
    gamma_ff_isdf/        % fdep=2,  isisdf=1  (+ local SAVE/)
    gamma_ff_dir/         % fdep=2,  isisdf=0
    gamma_cohsex_isdf/    % fdep=-2, isisdf=1
    gamma_cohsex_dir/     % fdep=-2, isisdf=0
```

Each case namelist uses:

- `groundstate_dir = '../qe.save'`
- `storage_dir = './SAVE'`

## Usage

```matlab
cd <repo>/examples
run_all
```

Or SLURM / batch:

```bash
cd examples
mkdir -p log
sbatch run_batch.sh
# or: bash run_batch.sh
```

`run_all` cleans each case’s local outputs (`SAVE/`, logs, `qp*.dat`), then runs
all four. QP tables land as `cases/<name>/qp*.dat`.
