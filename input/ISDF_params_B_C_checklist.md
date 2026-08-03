# ISDF 参数清单：B / C 类（待逐项检查）

对照规范：

- 变量名 → `input/allowed_param_list.m` 的 `def.ISDF`
- 默认值 → `input/default_param_values.m` 的 `def.ISDF`
- 调用处只读 `config.ISDF.<name>`，不在算法中间再写一套默认值

本文件只列 **B（orphan：未登记却硬编码/半 config）** 与 **C（已登记但不一致或几乎不用）**。  
A 类（已登记、但代码里还有 resolve fallback）不在此表。

检查时在 `[ ]` 改成决策，例如：

- `[x] 登记` — 写入 allowed + default，调用处去掉本地默认
- `[x] 保持硬编码` — 算法常量，不进 config
- `[x] 删除/废弃` — 参数或代码路径去掉
- `[x] 补 allowed` — 已有 default，补进 allowed

---

## B. Orphans（未在 allowed 登记，或根本不是 config 键）

### B1. 条件数 / Schur / 数值门槛

| # | 名称（建议键名） | 当前硬编码 | 位置 | 是否读 config | 决策 |
|---|------------------|------------|------|---------------|------|
| B1.1 | `max_cond_number`（Schur `init` 兜底） | `1e12` | `service/+isdf/+adaptive_double/isdf_schur_update.m`（及 single/adaptive 克隆） | 可选 varargin；正常由 `adaptive_max_cond_*` 传入 | [ ] |
| B1.2 | `condest` 门限 | `1e12` | `service/+isdf/cohsex_vcVnn.m` | 否 | [ ] |
| B1.3 | chol `jitter` | `~1e-12 * trace/N`（下限 `1e-12`） | `+adaptive_double/isdf_schur_update.m` | 否 | [ ] |
| B1.4 | `max_cch_cond` | `1e10` | `+coeff/gen_coeff_pseudo.m` | 否 | [ ] |
| B1.5 | `cnd` warn 门槛 | `1e6` | `+coeff/coeff_coarse_wf_filter.m` | 否 | [ ] |

### B2. 驱动 / 导入

| # | 名称（建议键名） | 当前硬编码 | 位置 | 是否读 config | 决策 |
|---|------------------|------------|------|---------------|------|
| B2.1 | backend 路由 `cutoff` | `1e-6` | `service/+isdf/driver.m`（single vs double） | 否 | [ ] |
| B2.2 | `import_isdf_adaptive_root` | `''` | `+adaptive_double/adaptiveisdf.m` / `+adaptive_single/adaptiveisdf.m` | **是**（`cfg_isdf.import_isdf_adaptive_root`），但 **未进 allowed/default** | [ ] |

### B3. 容差 / 对称 / 过滤

| # | 名称（建议键名） | 当前硬编码 | 位置 | 是否读 config | 决策 |
|---|------------------|------------|------|---------------|------|
| B3.1 | `rel_tol`（exclude） | `1e-5` | `+adaptive_*/isdf_exclude_point.m` | 否 | [ ] |
| B3.2 | `amp_threshold` | `1e-8` | `+coeff/coeff_coarse_wf_filter.m` | 否 | [ ] |
| B3.3 | `rot_tol` | `1e-4` | `+rsymm/bundle_refresh.m` | 否 | [ ] |
| B3.4 | kbz rot `tol` | `1e-5` | `+adaptive*/adaptiveisdf_kbz_rot_index.m` | 否（可回落到 lattice tol） | [ ] |
| B3.5 | `rgrid_symm_orbit` tol | `1e-3` | `rgrid_symm_orbit.m` | 否 | [ ] |
| B3.6 | `isdf_r_sampling_symm_orbit` tol | `1e-4` | `isdf_r_sampling_symm_orbit.m` | 否 | [ ] |

### B4. 占据数 / 校验报告

| # | 名称（建议键名） | 当前硬编码 | 位置 | 是否读 config | 决策 |
|---|------------------|------------|------|---------------|------|
| B4.1 | `tol_occ` / 占据划分 | `1e-6` 或 `1e-5` | `resolve_nrange.m` / `isdf_get_nrange.m` / `+validation/isdf_validate_HF*.m` | 否 | [ ] |
| B4.2 | HF mismatch 阈值 | `1e-4` | `+validation/isdf_validate_HF*.m`：默认一行摘要 (`rs`)；逐点仅 `v2l` | 否 | [ ] |

### B5. Gamma / SVD 路径上的“第二默认”（与已登记 `inv_*` 冲突）

> 若你决定以 `default_param_values` 为准，这些应归入“调用处禁止再写默认”，不一定是新键。

| # | 名称 | 当前硬编码 | 位置 | 与中央关系 | 决策 |
|---|------|------------|------|------------|------|
| B5.1 | `s_cut` 本地默认 | `0`（缺字段时） | `gen_tildeVq.m` / `prod_C_inv_*` | 中央已有 `inv_param` | [ ] |
| B5.2 | Gamma `inv_ratio` 本地 | 曾出现 `0.75` 等 | Gamma attach / `gen_tildeVq` 路径 | 中央已有 `inv_ratio` | [ ] |

---

## C. 已登记但不一致 / 几乎不用

### C1. 有 default、缺 allowed

| # | 键名 | 当前 default | 使用处 | 决策 |
|---|------|--------------|--------|------|
| C1.1 | `debug_checks` | `false` | `+debug/cache.m` / `on.m` | [ ] |
| C1.2 | `debug_level` | `'error'` | `+debug/cache.m` / `level.m` / `react.m` | [ ] |
| C1.3 | `debug_tags` | `[]` | `+debug/cache.m` / `on.m` | [ ] |

### C2. allowed + default 都有，但 `+isdf` 主路径几乎不读

| # | 键名 | 当前 default（以 `default_param_values.m` 为准） | 备注 | 决策 |
|---|------|--------------------------------------------------|------|------|
| C2.1 | `seed` | `0` | 旧 ISDF/Cauchy 风格；`+isdf` 未见读取 | [ ] |
| C2.2 | `init` | `'random'` | 同上 | [ ] |
| C2.3 | `weight` | `'add'` | 同上（勿与 `adaptive_weight` 混淆） | [ ] |
| C2.4 | `sys` | `[]` | 同上 | [ ] |
| C2.5 | `is_helper` | `false` | 同上 | [ ] |
| C2.6 | `inv_strategy` | `'dir'` | 登记了；`+isdf` 主路径是否消费待确认 | [ ] |
| C2.7 | `auto_inv_param` | `false` | 同上 | [ ] |
| C2.8 | `order` | `1` | 同上 | [ ] |
| C2.9 | `chol_maxit` | `10` | `iterated_cholqr` 有参默认 `10`，未必读 `config.ISDF.chol_maxit` | [ ] |

### C3. 已登记、有用，但方法字符串仍在代码里二次默认（偏 A/C 边界）

| # | 键名 | 当前 default | 代码二次默认 | 决策 |
|---|------|--------------|--------------|------|
| C3.1 | `exxmethod_type1` | `''` | vc → `'coarse'`（`+coeff/gen_coeff.m`） | [ ] |
| C3.2 | `exxmethod_type2` | `''` | vn → `'pseudo'` | [ ] |
| C3.3 | `exxmethod_type3` | `''` | nn → `'coarse'` | [ ] |
| C3.4 | `exxmethod` | `''` | 通用回落 | [ ] |

---

## 检查进度（可选总览）

- B 项合计：约 18（B1–B5）
- C 项合计：约 16（C1–C3）
- 完成：___ / ___

## 备注

1. `adaptive_max_cond_*` 等 A 类不在本表；你已在 `default_param_values.m` 改成 `1e12` 的，检查 Schur 兜底（B1.1）时要一起对齐。
2. `+adaptive` / `+adaptive_single` / `+adaptive_backends` 里常有同名硬编码克隆，决策时按“逻辑一处、实现多份”处理即可。
3. 本清单生成时对照的是 `service/+isdf` + 当前 `input/*_param_*.m`；若你改了 default，以仓库现状为准更新“当前 default”列。
