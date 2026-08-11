# util / common / input / interfaces 归类

日期：2026-08-11

---

## 当前库存

| 目录 | 文件 |
|------|------|
| `util/` | `constant_map`, `filename_map`, `allowed_param_list`, `default_param_values` |
| `common/` | `do_FFT`, `put_into_fftbox`, `generate_frequency`, `GaussLegendre`, `pseudoinv`, `testmemory` |
| `input/` | `read_input_param`, `validate_required_params`, `set_default_param_value`, `display_input_summary` |
| `interfaces/` | 见下（外部基态 → 内部 `data`） |

已迁出其它位置：

| 原位置 | 现位置 |
|--------|--------|
| `common/bz2ibz.m` | `service/+lattice/fold_k_to_BZ.m` |
| `common/find_gvec_in_glist.m` | `service/+lattice/find_gvec_in_glist.m` |
| `input/*ISDFCauchy*` | `isdf.Cauchy.*` |

---

## 职责归类

### A. 约定 / 表 → `util/`

maps + allow-list + defaults。

### B. 输入装配 → `input/`

namelist：parse → validate → defaults → summary。  
`isdf.Cauchy.setISDFCauchy` 由 `input_driver` 调用，不在 `input/`。

### C. 基态导入 → `interfaces/`（独立于 `service/`，对齐 Yambo）

```text
interfaces/
  load_groundstate_info.m      % C0 分发
  common/construct_rhoG.m      % C2
  qe/                          % C1（发布支持）
  kssolv/                      % C1（发布支持）
  formal/                      % 代码保留；分发器故意不开放
```

`QPstartup` 已 `add_mpaths_only(.../interfaces/)`。

### D. 数值小工具 → `common/`

含 **`generate_frequency` + `GaussLegendre`**（频率选点，可复用；不是 namelist 装配）。

---

更完整的用户向说明见 [`doc_rev/README.md`](../doc_rev/README.md) 与 [`doc_rev/class_reference.md`](../doc_rev/class_reference.md)。
