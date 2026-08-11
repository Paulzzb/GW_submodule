# ISDF report 命名与落盘规范

本文档整理 `+isdf` 下类 report 输出的命名约定与 `+report` 包结构。  
规则与全局约定一致：先在 `filename_map` 定义名字 → 在对应输出函数中读取 → 再输出。

---

## 1. 产物一览

| 产物 | 磁盘路径 | 写手 | 通道 |
|------|----------|------|------|
| HF 校验报告 | `isdf_report/o-ISDF_HF_id%d` | `isdf.report.hf` | `output.open` + `filename_map` |
| Adaptive phase-1 | `isdf_report/o-ISDF_adaptive_id%d` | `isdf.report.adaptive` | `output.open` + `filename_map` |
| 数值条件诊断 | `isdf_report/o-ISDF_cond` | `isdf.report.cond` | `output.open` + `filename_map` |
| ISDF run 摘要 | `r-<prefix>.log` | `run_summary` | `output.msg('r'/'rs')` |

System / Symmetry / FFT / Lattice / Wave functions 各段由对应 `*.driver` 结束时写入 `r-*`（不在 `+isdf/+report`）。

目录名固定为 `filename_map().isdf_report_dir`（`isdf_report`）；`isdf.driver` 入口若不存在则 `mkdir`。
| Coarse grid 摘要 | （无文件） | `+coeff/print_coarse_grid_report` | `fprintf` 屏幕 |
| Adaptive 填表 | （无文件） | `isdf.report.fill_adaptive` | 填 struct，供 adaptive writer 使用 |

全局 `r-*` / `l-*`（`output.msg('r'/'l')`）属于 run 级 report/log，**不是** ISDF 专属 OF。

旧名 `isdf_validate_HF_id*.txt` / `adaptiveisdf_id*.txt` / `numerical_cond_report` 入口已废弃并移除。

---

## 2. 统一规则（a / b / c）

1. **a.** 在 `util_profile/filename_map.m` 定义 basename 或模板。  
2. **b.** 在 `service/+isdf/+report/` 的输出函数中调用 `filename_map()` 取名。  
3. **c.** 一律经 `output.open` + `output.msg('o <logical_name>', …)` 落盘（不再手写 `fopen` 管 OF）。

边界：

| 内容 | 通道 |
|------|------|
| 专属诊断 / 校验产物 | 命名 OF：`o-ISDF_*` |
| 跑次参数 / 摘要 | 全局 report：`'r'` / `'rs'` |
| 迭代流水 / 调试 | `'l'` / `'v2l'` |

不要把 per-iq 条件数等大段明细塞进主 `r-*`。

---

## 3. 磁盘文件名

统一 Yambo 风格前缀 **`o-ISDF_<kind>`**。

| kind | `filename_map` 键 | 磁盘名 / 模板 | 实例策略 |
|------|-------------------|---------------|----------|
| 条件数 | `cond_report` | `o-ISDF_cond` | 一次 run 一份，可分段 append |
| HF 校验 | `hf_report` | `o-ISDF_HF_id%d` | 每个 pool id 一份 |
| Adaptive | `adaptive_report` | `o-ISDF_adaptive_id%d` | 每个 coarse id 一份 |

说明：

- **纵向/全程一份**（cond）：无 id 后缀。  
- **每个 pool 一份**（HF / adaptive）：保留 `_id%d`，便于多 slot；模板写在 `filename_map`，writer 内 `sprintf`。

---

## 4. `output.open` 逻辑名

逻辑名与磁盘 basename **相同**，一律取自 `filename_map` 的值（不再另写短名）：

| `filename_map` 键 | 磁盘 / OF 名 |
|-------------------|--------------|
| `cond_report` | `o-ISDF_cond` |
| `hf_report` | `o-ISDF_HF_id%d`（`sprintf` 后） |
| `adaptive_report` | `o-ISDF_adaptive_id%d`（`sprintf` 后） |

示例：

```matlab
def = filename_map();
of = def.cond_report;
output.open(of, fullfile(def.isdf_report_dir, of), 'w');
output.msg(['o ' of], '...');
```

带 id 的产物：

```matlab
def = filename_map();
of = sprintf(def.hf_report, id);
fpath = fullfile(def.isdf_report_dir, of);
output.open(of, fpath, 'w');
output.msg(['o ' of], '...');
```

---

## 5. `+report` 包结构

`+report` 只负责 **落盘专属 OF**；屏幕摘要走 `output.msg`，不必放进本包。

```
service/+isdf/+report/
  NAMING.md
  run_summary.m
  hf.m
  cond.m
  adaptive.m
  fill_adaptive.m
```

| 调用入口 | 说明 |
|----------|------|
| `isdf.report.run_summary` | 主 log ISDF run summary |
| `isdf.report.hf` | HF 校验报告 |
| `isdf.report.cond` | 数值条件诊断 |
| `isdf.report.adaptive` | Adaptive phase-1 报告 |
| `isdf.report.fill_adaptive` | 填 report struct（无文件） |
| `+coeff/print_coarse_grid_report` | **不在本包**；可改为 `output.msg('rs', …)` |

`+adaptive` 只保留算法；写报告只调 `isdf.report.*`。

---

## 6. `filename_map` 条目

```matlab
def.isdf_report_dir = 'isdf_report';
def.cond_report = 'o-ISDF_cond';
def.hf_report = 'o-ISDF_HF_id%d';
def.adaptive_report = 'o-ISDF_adaptive_id%d';
```

---

## 7. 迁移步骤（已完成）

1. ✅ 在 `filename_map` 补齐 `hf_report` / `adaptive_report`。  
2. ✅ 实现 `+report/cond.m`；调用点改为 `isdf.report.cond`。  
3. ✅ 改写 `+report/hf.m`；经 `output.open` 写 `o-ISDF_HF_id%d`。  
4. ✅ 迁入 `+report/adaptive` + `fill_adaptive`；`adaptiveisdf` 改调新入口。  
5. ✅ 测试 / smoke / `run_adaptiveisdf` 只认 `o-ISDF_*`。  
6. ✅ 删除 shim（`gen_report`、`numerical_cond_report`、`adaptiveisdf_write_phase1_report`、`adaptive_fill_report`）与旧名拷贝。

---

## 8. 已确认

1. HF / adaptive **保留** `_id%d`。  
2. 磁盘前缀一律 **`o-ISDF_`**。  
3. 函数用短名：`hf` / `cond` / `adaptive`。  
4. 兼容层（shim + 双写旧名）已移除。
