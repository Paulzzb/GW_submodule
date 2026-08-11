# Namelist 输入说明（用户版）

填写 case 目录下的 `test`（namelist）时用本文。只列跑通展示算例常用的参数。  
默认能量单位：**Ry**（`FREQUENCY.broadening` 除外，见下表）。

English version: [`../doc_user_EN/GW_input_description.md`](../doc_user_EN/GW_input_description.md)

ISDF 类型后缀：

| 后缀 | 对 | 含义 |
|------|----|------|
| `type1` | `vc` | 价带–导带 |
| `type2` | `vn` | 价带–价带 |
| `type3` | `nn` | 导带–导带 |

仓库示例：`examples/cases/*/test`。

---

## `&CONTROL`

| 参数 | 默认 | 说明 |
|------|------|------|
| `groundstate_dir` | （必填） | 基态目录，如 `'./qe.save'` |
| `groundstate_type` | `'qe'` | 当前暂只支持 `'qe'` |
| `storage_dir` | `'./SAVE/'` | 中间量目录；示例常用 `'./SAVE'` |
| `output_dir` | `'./'` | 报告与 `qp.dat` 所在目录 |
| `prefix` | `'QP'` | 报告文件前缀 → `r-<prefix>.log` |
| `log_level` | `1` | 日志详细程度 |
| `enable_k_points` | `0` | 多 k；Gamma 展示算例为 `.false.` / `0` |
| `isgw` | `1` | 启用 GW |
| `isbse` | `0` | BSE（当前不支持） |

---

## `&SYSTEM`

| 参数 | 默认 | 说明 |
|------|------|------|
| `number_bands_in_summation` | `-1` | 求和用能带数；`-1` → 由基态填充 |
| `energy_band_index_min` | `-1` | QP 输出能带下界 |
| `energy_band_index_max` | `-1` | QP 输出能带上界 |

---

## `&CUTOFFS`

| 参数 | 默认 | 说明 |
|------|------|------|
| `coulomb_truncation_method` | `2` | `0` 无截断；`2` 球截断（常用） |
| `coulomb_truncation_parameter` | `5.0` | 截断参数（Ry；method `2` 时为球半径） |
| `coulomb_cutoff` | `-1.0` | Coulomb/G 截断（Ry）；`-1` → 用基态波函数截断 |
| `density_cutoff` | `-1` | 密度截断（Ry）；`-1` → 约 2× 波函数截断；当前无用 |

截断取值敏感，改动前请对照基态与示例。

---

## `&FREQUENCY`

| 参数 | 默认 | 说明 |
|------|------|------|
| `frequency_dependence` | `0` | **`-2`** COHSEX（Gamma）；**`2`** 全频 CD（Gamma）。展示算例主要用这两档 |
| `frequency_dependence_method` | `2` | 全频时的方法选择 |
| `frequency_low_cutoff` | `-1.0` | 低频截止（Ry）；`-1` → 由能带估计 |
| `broadening` | `0.25` | 全频展宽（**eV**） |
| `delta_frequency` | `≈0.147` | 实频网格步长（Ry） |
| `number_imaginary_freqs` | `15` | 虚频点数 |
| `cd_integration_parameter` | `2.0` | CD 积分参数（Ry） |

`qp.launcher` 路由：

| `frequency_dependence` | 计算路径 |
|------------------------|----------|
| `-2` | `gw.x_Gamma` + `gw.cohsex_Gamma` |
| `2` | `gw.fullfreq_cd_res_Gamma` + `gw.fullfreq_cd_int_Gamma` |

---

## `&ISDF`

| 参数 | 默认 | 说明 |
|------|------|------|
| `isisdf` | `1` | `1` 用 ISDF；`0` 走 dense（对比算例） |
| `compute_vc` / `compute_vn` / `compute_nn` | `true` | 是否构建对应类型的 ISDF |
| `isdf_ratio_type1/2/3` | `8.0` | 各类型 ISDF 比率（须 &gt; 0） |
| `exxmethod` | → `'pseudo'` | 插值点方法：`'pseudo'` / `'qrcp'` / `'kmeans'` / `'coarse'` 等 |
| `exxmethod_type1/2/3` | 空 → 跟 `exxmethod` | 可按类型覆盖 |
| `iscauchy` | `false` | 是否启用 Cauchy 积分相关选项 |
| `validate_hf` | `false` | 是否写 HF 类校验报告（`isdf_report/`） |
| `is_helper` | `false` | 辅助量 / 旁路 |

自适应细化（`adaptive_*`）有默认值，一般沿用示例即可；细调时对照 `examples/cases/*/test`。

目标秩规模大致：\(N_\mu \approx \texttt{isdf\_ratio} \sqrt{N_1 N_2}\)。

---

## `&COHSEX`

| 参数 | 默认 | 说明 |
|------|------|------|
| `exact_ch` | `false` | 是否用精确 COH 求和（仅支持 `isisdf=1`） |
| `ex_use_which_isdf` | `'vn'` | 交换用哪路 ISDF：`'vn'` 或 `'nn'` |

在 `frequency_dependence = -2` 时本块才主要相关。

---

## 最小示例（QE + 全频 ISDF）

```text
&CONTROL
  groundstate_dir = '../qe.save',
  groundstate_type = 'qe',
  storage_dir = './SAVE',
  output_dir = './',
  prefix = 'ff_isdf',
  log_level = 1,
END &CONTROL

&SYSTEM
  number_bands_in_summation = 319,
  energy_band_index_min = 1,
  energy_band_index_max = 32,
END &SYSTEM

&CUTOFFS
  coulomb_truncation_method = 2,
  coulomb_truncation_parameter = 5.0,
  coulomb_cutoff = 15.0,
END &CUTOFFS

&FREQUENCY
  frequency_dependence = 2,
  frequency_dependence_method = 2,
  frequency_low_cutoff = 25.0,
  broadening = 0.01876,
  number_imaginary_freqs = 25,
END &FREQUENCY

&ISDF
  isisdf = 1,
  isdf_ratio_type1 = 12.0,
  isdf_ratio_type2 = 24.0,
  isdf_ratio_type3 = 16.0,
  exxmethod = 'qrcp',
END &ISDF
```

更多完整示例见 `examples/cases/`。
