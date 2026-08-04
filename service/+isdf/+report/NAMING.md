# ISDF report 命名与落盘规范

本文档整理 `+isdf` 下类 report 输出的现状、统一命名约定，以及迁入 `+report` 的目标结构。  
规则与全局约定一致：先在 `filename_map` 定义名字 → 在对应输出函数中读取 → 再输出。

---

## 1. 现状盘点

| 产物 | 磁盘名（现状） | 写手 | 通道 | 备注 |
|------|----------------|------|------|------|
| HF 校验报告 | `o-ISDF_HF_id%d`（兼容写旧名） | `+report/hf`（`gen_report` 为 shim） | `output.open` + `filename_map` | 步骤 3 已完成 |
| Adaptive phase-1 | `adaptiveisdf_id<coarse_id>.txt` | `+adaptive/adaptiveisdf_write_phase1_report` | 手写 `fopen` | single/double 共用 |
| 数值条件诊断 | `o-ISDF_cond` | `+report/cond`（旧入口 `numerical_cond_report` 为 shim） | `output.open` + `filename_map` | 步骤 2 已完成 |
| Coarse grid 摘要 | （无文件） | `+coeff/print_coarse_grid_report` | `fprintf` 屏幕 | 名字带 report，但不是文件产物 |
| Adaptive 填表 | （无文件） | `+adaptive/adaptive_fill_report` | 填 struct | 供 phase-1 writer 使用 |

全局 `r-*` / `l-*`（`output.msg('r'/'l')`）属于 run 级 report/log，**不是** ISDF 专属 OF。

问题：三套真正的文件产物、三套命名风格、两套 IO 路径；`+report` 目前只有 HF 的 `gen_report`。

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

统一 Yambo 风格前缀 **`o-ISDF_<kind>`**，与已有 `o-ISDF_cond` 对齐。

| kind | 建议 `filename_map` 键 | 建议磁盘名 / 模板 | 实例策略 |
|------|------------------------|-------------------|----------|
| 条件数 | `cond_report`（已有） | `o-ISDF_cond` | 一次 run 一份，可分段 append |
| HF 校验 | `hf_report` | `o-ISDF_HF_id%d` | 每个 pool id 一份 |
| Adaptive | `adaptive_report` | `o-ISDF_adaptive_id%d` | 每个 coarse id 一份 |

说明：

- **纵向/全程一份**（cond）：无 id 后缀。  
- **每个 pool 一份**（HF / adaptive）：保留 `_id%d`，便于多 slot；模板写在 `filename_map`，writer 内 `sprintf`。  
- 旧名 `isdf_validate_HF_id*.txt` / `adaptiveisdf_id*.txt` 废弃；测试 glob 需同步更新。

---

## 4. `output.open` 逻辑名

逻辑名与磁盘 basename 分开，用稳定短名：

| 逻辑名（`how`） | 磁盘 |
|-----------------|------|
| `cond_report` | `o-ISDF_cond` |
| `hf_report` | `o-ISDF_HF_id%d` |
| `adaptive_report` | `o-ISDF_adaptive_id%d` |

示例：

```matlab
def = filename_map();
output.open('cond_report', fullfile(pwd, def.cond_report), 'w');
output.msg('o cond_report', '...');
```

带 id 的产物：

```matlab
def = filename_map();
fpath = fullfile(outDir, sprintf(def.hf_report, id));
output.open('hf_report', fpath, 'w');
output.msg('o hf_report', '...');
```

---

## 5. `+report` 包结构与函数命名

`+report` 只负责 **落盘专属 OF**；屏幕摘要走 `output.msg`，不必放进本包。

目标结构：

```
service/+isdf/+report/
  NAMING.md              % 本文档
  hf.m                   % 原 gen_report
  cond.m                 % 原 numerical_cond_report
  adaptive.m             % 原 adaptiveisdf_write_phase1_report
  fill_adaptive.m        % 原 adaptive_fill_report（可选，纯填表）
```

| 现状 | 迁入后 |
|------|--------|
| `+report/gen_report` | `isdf.report.hf` |
| `numerical_cond_report`（`+isdf` 根） | `isdf.report.cond` |
| `+adaptive/adaptiveisdf_write_phase1_report` | `isdf.report.adaptive` |
| `+adaptive/adaptive_fill_report` | `isdf.report.fill_adaptive`（可选） |
| `+coeff/print_coarse_grid_report` | **不迁入**；改为 `output.msg('rs', …)` |

`+adaptive` 只保留算法；写报告只调 `isdf.report.*`。

函数名偏好（已倾向短名）：

- `report.hf` / `report.cond` / `report.adaptive`  
- 备选更明示：`write_hf` / `write_cond` / `write_adaptive`

---

## 6. `filename_map` 建议条目

```matlab
def.cond_report = 'o-ISDF_cond';
def.hf_report = 'o-ISDF_HF_id%d';
def.adaptive_report = 'o-ISDF_adaptive_id%d';
```

（三项均已写入 `filename_map.m`。）

---

## 7. 迁移步骤（建议顺序）

1. ✅ 在 `filename_map` 补齐 `hf_report` / `adaptive_report`。  
2. ✅ 实现 `+report/cond.m`；调用点改为 `isdf.report.cond`；`isdf.numerical_cond_report` 保留为兼容 shim。  
3. ✅ 改写 `+report/hf.m`；`gen_report` 为 shim；兼容拷贝旧名 `isdf_validate_HF_id*.txt`；测试 glob 认新旧两种。  
4. 迁入 adaptive phase-1 writer + fill；`adaptiveisdf` 改为调 `isdf.report.adaptive`；测试一轮兼容旧名 `adaptiveisdf_id*.txt`。  
5. 更新测试 / smoke 中的旧 glob，去掉兼容层。  
6. 删除旧路径文件与过时注释。

每步独立验收后再动下一步。

---

## 8. 已确认

1. HF / adaptive **保留** `_id%d`。  
2. 磁盘前缀一律 **`o-ISDF_`**。  
3. 函数用短名：`hf` / `cond` / `adaptive`。  
4. 迁移期 **先一轮兼容**（shim + 测试仍可认旧文件名），再清旧名。
