# examples/ — public ISDF vs dense demo

Small Gamma-Si showcase: same groundstate, compare **ISDF** vs **dense** for
full-frequency (`fdep=2`) and COHSEX (`fdep=-2`).

Internal verification suites live under `tests/` (not this tree).

## Layout

```text
examples/
  run_all.m / run_all_profiled.m / run_all_shared_save.sh
  run_cohsex.m
  cases/
    qe.save/              % shared QE groundstate
    SAVE/                 % shared storage (created at runtime)
    gamma_ff_isdf/        % fdep=2,  isisdf=1
    gamma_ff_dir/         % fdep=2,  isisdf=0
    gamma_cohsex_isdf/    % fdep=-2, isisdf=1
    gamma_cohsex_dir/     % fdep=-2, isisdf=0
```

Each case namelist uses:

- `groundstate_dir = '../qe.save'`
- `storage_dir = '../SAVE'`

## Usage

```matlab
cd <repo>/examples
run_all
collect_qp_energies
```

Or SLURM:

```bash
cd examples
mkdir -p log
sbatch run_all_shared_save.sh
```

`run_all` clears shared `cases/SAVE` once, then runs ISDF cases before dense.
