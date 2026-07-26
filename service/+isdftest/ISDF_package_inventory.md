# `+isdf` 包函数与职责一览

> 生成目的：为后续**重命名**与**分类/子目录拆分**提供清单。  
> 路径：`GW/service/+isdf/`（子包：`+isdf/+base/`）。  
> 对外调用形式：**根目录**为 **`isdftest.<文件名>`**；**子包**为 **`isdftest.<子包>.<文件名>`**（见下）。

---

## 0. 子包布局（2026-04 起）

| 子目录 | MATLAB 前缀 | 内容 |
|--------|----------------|------|
| `+isdf/+coeff/` | `isdftest.coeff.*` | 产生/使用 `coeff_seper` 的粗网格管线：`gen_coeff*`、`coeff_coarse_*`、`isdf_coarse_R_rot`、`divisors_int32`、`print_coarse_grid_report`、`isdf_get_coeff`、`isdf_apply_symm_on_coarse` 等 |
| `+isdf/+adaptive/` | `isdftest.adaptive.*` | 自适应增心：`adaptiveisdf`、`adaptive_weight`、`isdf_schur_update`、`isdf_prod`、`run_adaptiveisdf`、`isdf_exclude_point`、`adaptiveisdf_kbz_rot_index` |
| `+isdf/+validation/` | `isdftest.validation.*` | 效果验证：`isdf_validation`、`isdf_validate_HF` |
| `+isdf/+report/` | `isdftest.report.*` | 文本报告：`gen_report` |
| `+isdf/` 根目录 | `isdftest.*` | 池与门面：`manager`、`get`、`save2mod`、`free`、`driver`、`isdftest_add`、`gen_tildeVq`、`isdftest_get_nrange`、R 网格轨道辅助等（其余待你继续整理） |
| `+isdf/+base/` | `isdftest.base.*` | 数据类 `isdf_m` |

---

## 1. 架构速览

- **对象池**：最多 `N_MAX=10` 个 `isdftest.base.isdftest_m` 实例，由 `manager` 持久化管理；`get` / `save2mod` / `free` / `isdftest_add` 等是薄封装。
- **粗网格 ISDF 主路径**：`isdftest.coeff.gen_coeff` → `gen_coeff_coarse` → `coeff_coarse_fft_grid_build` + `coeff_coarse_wf_extract` → 写回 `coeff_seper`；再 `isdftest.gen_tildeVq`。
- **校验**：`isdftest.validation.isdf_validation` 按 `interp_scheme` 分发；HF 细节在 `isdftest.validation.isdf_validate_HF` + `isdftest.report.gen_report`。
- **自适应**：`isdftest.adaptive.adaptiveisdf` 依赖同包内 `isdf_schur_update` / `adaptive_weight` / `isdf_prod` 及根目录 `isdftest.isdftest_get_nrange` 等。

---

## 2. 按职责分类

### 2.1 对象池与当前指针

| 函数文件 | 功能 |
|----------|------|
| `manager` | 核心：池的 `get` / `save2mod` / `add` / `select` / `current` / `nmax` / `list` / `free`；校验 id 的局部函数 `isdftest_manager_validate_id`。 |
| `get` | 读指定 id 或当前 id 的 `isdf_m`，不改变 current（除非省略 id 时用 current）。 |
| `save2mod` | 写回池；可指定 id。 |
| `free` | 清空整个池或单个 id。 |
| `isdftest_add` | 分配第一个空槽，设 `desc`，设为 current，返回 id。 |
| `isdftest_list` | 返回各槽的 `id` / `desc` / `assigned` / `empty` 结构体数组。 |
| `isdftest_current` | 当前池 id（int32）。 |
| `isdftest_nmax` | 池容量 `N_MAX`。 |
| `isdftest_select` | 设置 current id（转调 `manager('select', …)`）。 |

### 2.2 与 relay / 测试入口相关的 driver

| 函数文件 | 功能 |
|----------|------|
| `driver` | **GW relay 用**：按 `config.ISDF` 建 vn/nn 槽、`gen_coeff`、`print_coarse_grid_report`、`gen_tildeVq`、（可选）校验。 |
| `isdftest_driver` | **独立 input 目录**：读 `GWinput`/`config`，写 `ISDF_DB`；旧式一站式 ISDF 驱动。 |
| `run_adaptiveisdf` | 在**当前 profile**（需 `SAVE/config.mat`、`test_relay_stage.mat`）下重建 coarse 槽、`gen_coeff`，再 `cd` 到 `output_dir` 跑 `adaptiveisdf` 并收集报告路径。 |

### 2.3 系数与粗网格管线（`interp_scheme == "coarse"`）

| 函数文件 | 功能 |
|----------|------|
| `gen_coeff` | **路由**：据 `config.exxmethod` 或 `isdfoptions.samp` 选择 `default` / `qrcp` / `kmeans` / `coarse`；局部函数 `local_pick_method`。 |
| `gen_coeff_coarse` | **当前主实现**：`coeff_coarse_fft_grid_build` → 取 `R_sampling_RLU`/`R_rot_coarse` → `coeff_coarse_wf_extract` → 写 `coeff_seper`、`interp_scheme="coarse"` 等。 |
| `coeff_coarse_fft_grid_build` | 由 `nisdf` 与细网格推正则粗 FFT 盒 `fftgrid_c`，并准备粗网格相关量。 |
| `coeff_coarse_wf_extract` | 在粗 RLU 点上抽取波函数系数（FFT / G 表 / 映射）。 |
| `isdf_coarse_R_rot` | 在规则粗 FFT RLU 格点上构造对称性指标映射 `R_rot_coarse_*`（`nr × nsym`）。 |
| `gen_coarse_Rgrid` | 细网格上按 gcd/step 规则选粗子格点；`Nmu_tmp = isdf_ratio * nrep`；依赖 `divisors_int32`。 |
| `gen_indices_coarse` | **遗留名**：注释称与 `gen_coeff_coarse` 的 r-grid模式相同；实现为调用 `isdftest.gen_coeff_coarse()`（**无 id 参数，可能与当前 `gen_coeff_coarse(id)` 签名不一致**，重构时需核对）。 |
| `print_coarse_grid_report` | 打印密/粗 FFT 尺寸与点数比；非 coarse 则直接返回。 |
| `divisors_int32` | 正整数 `n` 的所有正因数（int32 列），供粗网格步长搜索。 |

### 2.4 占位 / 未实现的索引算法

| 函数文件 | 功能 |
|----------|------|
| `gen_coeff_default` | 占位：警告并返回空索引。 |
| `gen_coeff_qrcp` | 占位：QRCP 路径。 |
| `gen_coeff_kmeans` | 占位：k-means 路径。 |

### 2.5 q 空间与 EXX 型张量

| 函数文件 | 功能 |
|----------|------|
| `gen_tildeVq` | 由粗 ISDF 与 Coulomb 构建 `tildeVq`（及 helper 量）；内部需 `isdftest_get_nrange`。 |

### 2.6 散射/接口用的系数与对称性

| 函数文件 | 功能 |
|----------|------|
| `isdf_get_coeff` | 按 `interp_scheme` 返回 `(c1,c2)` 型系数或做 coarse 下对称性旋转；支持 `'reset'` 清持久量。 |
| `isdf_apply_symm_on_coarse` | 在粗格点上应用对称 `isym`（置换 ± 时间反演）；支持 `'reset'`。 |

### 2.7 校验与报告文件

| 函数文件 | 功能 |
|----------|------|
| `isdf_validation` | 根据 `interp_scheme`（coarse / adaptive）调用相应校验路径（与 `isdf_validate_HF` 配合）。 |
| `isdf_validate_HF` | **HF 能量对比**：直接交换 vs ISDF `tildeVq` 收缩；写结构化 `report`；可写 `isdf_validate_HF_id<id>.txt`；含进度与计时。 |
| `gen_report` | 将 `isdf_validate_HF` 的 `report` 写成文本文件；局部函数 `isdf_gen_report_table2char`。 |

**同文件内非导出辅助函数（`isdf_validate_HF.m`）**

- `isdf_validate_HF_stream_init` / `isdf_validate_HF_stream_push` / `isdf_validate_HF_stream_finalize`：流式统计。
- `isdftest_ensure_timing_initialized`：测试/计时初始化。

### 2.8 自适应 ISDF（增广插值点）

| 函数文件 | 功能 |
|----------|------|
| `adaptiveisdf` | 对给定 coarse **池 id** 做自适应加心、损失迭代、轨道报告、`isdf_schur_update` 校验；新建 `interp_scheme='adaptive'` 的池并 `gen_tildeVq` + `isdf_validation`。 |
| `adaptive_weight` | 与 Schur 更新配套的**权重**持久状态：`init` / `init_update` / `update` / `get` / `clear`。 |
| `isdf_schur_update` | **Gram / Cholesky 型**分块更新的持久状态：`init` / `update` / `get` / `clear`。 |
| `isdf_exclude_point` | 可选：按对角近似权重剔除病态 centroid（调试/预处理）。 |
| `isdftest_get_nrange` | 按池 `desc`（`nn`/`vn`/`vc`）返回 ISDF 用的能带索引范围 `nrange1`/`nrange2`。 |
| `isdf_prod` | 四组波函数块上的 `(Psi*psi').*(Phi*phi')` 型乘积（用于 Gram/权重）。 |

**同文件内局部函数（`adaptiveisdftest.m`）**

- `adaptiveisdf_write_phase1_report`：写 `adaptiveisdf_id<coarse_id>.txt`。
- `adaptiveisdf_orbit_mask_bfs` / `adaptiveisdf_orbit_union_bfs`：FFT 格点上对称轨道 BFS。
- `adaptiveisdf_r_sampling_rows_to_lin`：将 `R_sampling_RLU` 行映射到细格线指标。

### 2.9 实空间 R 网格与对称轨道（对比 / 诊断）

| 函数文件 | 功能 |
|----------|------|
| `rgrid_symm_orbit` | 在**全 FFT R 网格**上构造对称轨道分解（`ir2rep` / `ir2rot` 等）。 |
| `isdftest_r_sampling_symm_orbit` | 在**当前 ISDF 的 `R_sampling_RLU`** 上构造同类轨道结构；可选 id，默认 current。 |

**同文件内局部函数（`isdftest_r_sampling_symm_orbit.m`）**

- `local_mod_rows` / `local_periodic_dist` / `local_key`：周期边界与哈希键。

### 2.10 k 空间（全 BZ）

| 函数文件 | 功能 |
|----------|------|
| `adaptiveisdf_kbz_rot_index` | 全 BZ k 指标在晶体对称 `irot` 下的像（RLU，模倒格矢）；用于 k 展开一致性。 |

### 2.11 数据模型（类）

| 文件 | 功能 |
|------|------|
| `+base/isdf_m` | `classdef`：池内一条 ISDF 记录（`desc`、`nisdf`、`coeff_seper`、`tildeVq`、`helperqG`、`fftgrid_c`、`R_rot_coarse`、`R_rot_extra`、`interp_scheme`、`N_coarse`、`N_extra` 等）；构造函数默认 `allocated/assigned`。 |

---

## 3. 文件统计（约）

| 类型 | 数量（约） |
|------|------------|
| `+isdf/*.m` 顶层函数文件 | 40+ |
| `+isdf/+base/*.m` | 1（类） |
| 单文件内多个 `function`（局部/辅助） | `adaptiveisdftest.m`、`isdf_validate_HF.m`、`isdftest_r_sampling_symm_orbit.m`、`gen_report.m`、`gen_coeff.m`、`manager.m` 等 |

---

## 4. 命名与分类上的明显痛点（供你决断）

1. **`isdf_*` 前缀重复**：包名已是 `isdf`，函数名再 `isdftest_get_nrange`、`isdf_validation` 等冗长。
2. **`gen_*` 混用**：有的表示路由（`gen_coeff`），有的表示具体算法（`gen_tildeVq`），有的是占位（`gen_coeff_kmeans`）。
3. **两套 driver**：`driver`（relay）、`isdftest_driver`（旧目录流）、`run_adaptiveisdf`（自适应专用）职责交叉但入口不同。
4. **`isdf_validation` vs `isdf_validate_HF`**：名称接近，层级一个是总入口、一个是 HF 细节。
5. **`gen_indices_coarse` vs `gen_coeff_coarse`**：遗留别名，且与当前 `gen_coeff_coarse(id)` 签名需再核实。

---

## 5. 修订记录

| 日期 | 说明 |
|------|------|
| 2026-04-19 | 初版：按当前仓库 `+isdf` 源码整理。 |
| 2026-04-19 | 子包拆分：`+coeff` / `+adaptive` / `+validation` / `+report` 与调用名更新。 |
