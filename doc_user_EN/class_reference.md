# Result structures (user guide)

Users normally only need `E` returned by `qp.launcher` / `qp_driver`, and the text table `qp.dat`.  
Column meanings and file layout: [`output_description.md`](output_description.md).

中文版：[`../doc_user_ZH/class_reference.md`](../doc_user_ZH/class_reference.md)

---

## `E` (return value of `qp.launcher`)

| Field | Description |
|-------|-------------|
| `Eqp` | Quasiparticle energy: `Ex + Esx_x + Ecoh` (**Ry**, selected bands × k) |
| `Ex` | Exact exchange Σ_x (Ry) |
| `Esx_x` | Screened-exchange correction (COHSEX SX−X, or full-freq residual) (Ry) |
| `Ecoh` | Coulomb-hole / integral contribution (Ry) |
| `Eqp0` | Not filled by the launcher; table uses `Eo + Sig − Vxc` |
| `fout` | Path to written `qp.dat` (empty on failure) |

```matlab
E = qp.launcher(config);
% or
E = qp_driver('./SAVE');
```

In-memory arrays are **Ry**; values in `qp.dat` are **eV** (`ry2ev ≈ 13.6057`).

---

## `qp.dat`

- Path: usually under `CONTROL.output_dir`; file name fixed as `qp.dat`
- Written by `qp.fout` (called inside `qp.launcher`)
- Columns include band index and Emf / Eo / X / SX−X / CH / Sig / Vxc / Eqp0, …

Full headers and units: [`output_description.md`](output_description.md#4-qpdat).

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-11 | ZZ | Split into `doc_user_ZH` / `doc_user_EN` |
| 2026-08-11 | ZZ | User edition: keep only `E` / `qp.dat` |
