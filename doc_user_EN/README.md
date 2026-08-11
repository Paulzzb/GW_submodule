# Quick-start guide

This guide is for **users**: prepare a ground state, fill a namelist, run a calculation, and inspect results.

中文版：[`../doc_user_ZH/README.md`](../doc_user_ZH/README.md)

| Document | Content |
|----------|---------|
| [`GW_input_description.md`](GW_input_description.md) | Namelist parameters (user slim) |
| [`output_description.md`](output_description.md) | Output files |
| [`class_reference.md`](class_reference.md) | Result struct `E` / `qp.dat` |
| [`../examples/README.md`](../examples/README.md) | Small demos (ISDF vs dense) |

---

## 1. Directories you will use

```text
QP_root/
├── QPstartup.m         # Startup: paths + optional MEX build
├── driver/             # input_driver / qp_driver
├── examples/           # Showcase cases (start here)
├── doc_user_EN/        # This directory: English user docs
├── doc_user_ZH/        # Chinese user docs
├── interfaces/         # Load QE / KSSOLV ground state
├── input/              # Parse namelist
├── service/ · packages/ · util/ · common/ · …
└── …
```

Main runtime products:

- `SAVE/` (or your `storage_dir`): `data.mat`, `config.mat`, `relay_stage.mat`
- `qp.dat`, `r-<prefix>.log`
- Optional: `isdf_report/`

---

## 2. Workflow (sketch)

```text
  QE qe.save  or  KSSOLV groundstate.mat
              |
              v
         namelist (./test)
              |
              v
         input_driver
              |
      +-------+--------+
      v                v
  SAVE/data.mat   SAVE/config.mat
      |                |
      +-------+--------+
              v
         build services → relay_stage.mat
              |
              v
         qp.launcher / qp_driver
              |
              v
         E  +  qp.dat  +  r-*.log
```

`qp.launcher` selects the path from `FREQUENCY.frequency_dependence`:

| Value | Path |
|---:|------|
| `-2` | COHSEX (Gamma): `gw.x_Gamma` + `gw.cohsex_Gamma` |
| `2` | Full-frequency CD (Gamma): `gw.fullfreq_cd_res_Gamma` + `gw.fullfreq_cd_int_Gamma` |

---

## 3. Quick start

### Step 1: Initialize the environment

In MATLAB, from the repository root:

```matlab
QPstartup;
```

This resets the path and builds the ISDF MEX when needed.

### Step 2: Prepare the ground state (QE recommended)

- QE must support **HDF5** wavefunction output.
- Run `pw.x` `scf` (and `bands` / `nscf` if needed).
- Export `vxc.dat` with `pw2bgw.x` (or an equivalent workflow).
- Place the following in one directory (e.g. `./qe.save/` or `examples/cases/qe.save/`):

```text
charge-density.hdf5
data-file-schema.xml
wfc*.hdf5
vxc.dat
```

In the namelist:

```text
&CONTROL
  groundstate_dir = './qe.save',   % or '../qe.save'
  groundstate_type = 'qe',
  storage_dir = './SAVE',
  output_dir = './',
  prefix = 'Si_gamma',
END &CONTROL
```

**KSSOLV:** (not supported yet) place `groundstate.mat` (variable name `groundstate`) and set `groundstate_type = 'kssolv'`.

Parameters: [`GW_input_description.md`](GW_input_description.md). Runnable examples: `examples/cases/*/test`.

### Step 3: Run

Recommended (same as `examples/run_cohsex.m`):

```matlab
input_driver('./test');                 % build SAVE/ and services
load('./SAVE/config.mat', 'config');
E = qp.launcher(config);                % QP calculation
```

Two-stage (resume from stage in another MATLAB session):

```matlab
input_driver('./test');
E = qp_driver('./SAVE');
```

Showcase only:

```matlab
cd examples
run_all
```

### Step 4: Check results

| File | Meaning |
|------|---------|
| `r-<prefix>.log` | Run report |
| `qp.dat` | Quasiparticle table (eV) |
| `SAVE/data.mat` | Ground-state data |
| `SAVE/config.mat` | Parsed configuration |
| `SAVE/relay_stage.mat` | Service cache |
| `isdf_report/` | ISDF diagnostics (if enabled) |

See [`output_description.md`](output_description.md); fields of `E` in [`class_reference.md`](class_reference.md).

### Stage cache (brief)

- If `relay_stage.mat` and `data.mat` already exist, `input_driver` reuses `data` but **always** rebuilds `config` from the namelist.
- Changing CUTOFFS / ISDF usually does not require deleting the stage; if the ground state or symmetry-related settings change, delete `data.mat` and `relay_stage.mat` and rerun.

---

## 4. Related docs

| Document | Content |
|----------|---------|
| [`GW_input_description.md`](GW_input_description.md) | Namelist (user slim) |
| [`output_description.md`](output_description.md) | Output files |
| [`class_reference.md`](class_reference.md) | `E` / `qp.dat` |
| [`../examples/README.md`](../examples/README.md) | Showcase cases |

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-11 | ZZ | Split into `doc_user_ZH` / `doc_user_EN` |
| 2026-08-11 | ZZ | User surface: drop extending section; slim namelist / class |
| 2026-08-08 | Zhengbang | Rewrite for service/packages architecture |
