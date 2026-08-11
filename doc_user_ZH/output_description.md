# 输出说明

本文说明当前 QP 框架的主要输出文件。

相关文档：[`README.md`](README.md)、[`class_reference.md`](class_reference.md)。  
文件名常量见 `util/filename_map.m`。

English version: [`../doc_user_EN/output_description.md`](../doc_user_EN/output_description.md)

---

## 总览

| 位置 | 文件 | 写出方 |
|------|------|--------|
| `CONTROL.storage_dir`（如 `./SAVE`） | `data.mat`、`config.mat`、`relay_stage.mat`、ISDF checkpoint | `input_driver` / `service_driver` / ISDF |
| `CONTROL.output_dir`（如 `./`） | `r-<prefix>.log`、`qp.dat` | `display_input_summary` / `qp.fout` |
| case 工作目录 | `isdf_report/o-ISDF_*` | ISDF 报告 |

默认 `storage_dir` 为 `./SAVE/`；示例通常显式设为 `./SAVE`。

---

## 1. `data.mat`

**路径：** `<storage_dir>/data.mat`  
**写出：** `input_driver`  
**变量：** `data`

**说明：**  
经 `interfaces/load_groundstate_info` 载入的基态结构体（通常为 QE `*.save`）。

**用途：**  
持久化 DFT 输入，后续阶段不必再读 QE HDF5。

**格式：** MATLAB `-v7.3`（无压缩）。

---

## 2. `config.mat`

**路径：** `<storage_dir>/config.mat`  
**写出：** `input_driver`  
**变量：** `config`

**说明：**  
解析后的 namelist 加默认值（`set_default_param_value`）。  
常见块：`CONTROL`、`SYSTEM`、`CUTOFFS`、`FREQUENCY`、`ISDF`、`COHSEX`。

**重要：**  
`input_driver` **始终** 按 namelist 重建 `config`。  
即使保留 stage 缓存，修改 CUTOFFS / ISDF / FREQUENCY 后下次调用 `input_driver` 也会生效。

参数概览：[`GW_input_description.md`](GW_input_description.md)。

---

## 3. `relay_stage.mat`

**路径：** `<storage_dir>/relay_stage.mat`  
**写出：** `service_driver` 内的 `relay.save2db`  
**变量：** `relay_stage`

**说明：**  
缓存的 service 层对象（昂贵、由 data 导出）：

- 可恢复：`system`、`symmetry`、`FFT`、`wave_functions`（及相关）
- **不**进 stage（每次按 config 重建）：`lattice`、`coulomb`、ISDF

**用途：**  
让 `qp_driver` 或另一次 MATLAB 会话通过 `relay.stage_from_db` + `relay.restore` 跳过从头重建波函数。

若存在 `relay_stage.mat` 但缺少 `data.mat`，`input_driver` 会警告并从基态目录重建。

---

## 4. `qp.dat`

**路径：** `<output_dir>/qp.dat`（若未设 `output_dir` 则回退到 `storage_dir`）  
**写出：** `qp.fout`（由 `qp.launcher` 调用）  
**结构：** 结果结构体 `E` — 见 [`class_reference.md`](class_reference.md)

**说明：**  
可读的准粒子能带表。

**单位：** 表中能量均为 **eV**。

**注：**  
`CONTROL.outfile`（默认 `'GWoutput'`）目前 **不被** `qp.fout` 使用；文件名固定为 `qp.dat`。

### 静态 / COHSEX 表头

当 `frequency_dependence` 不是 `2`（例如 `-2`）时：

```plaintext
   n         Emf          Eo           X        SX-X          CH         Sig         Vxc        Eqp0
```

每个能带一行（k 循环顺序与 `qp.fout` 一致）。

### 全频表头（`frequency_dependence == 2`）

每个能带两行：第一行实部，第二行为 SX-X / CH / Sig / Eqp0 的虚部：

```plaintext
   n         Emf          Eo           X      Re SX-X       Re CH      Re Sig        Vxc     Re Eqp0
                                     Im SX-X       Im CH      Im Sig                Im Eqp0
  13    3.467416    3.467416  -12.445148    1.211312   -0.172058  -11.405895  -10.521226    2.582748 
                                            0.000059    0.000000    0.000059                0.000059 
```

### 列含义

| 列 | 含义 |
|----|------|
| `n` | 能带指标（`energy_band_index_min` … `max`） |
| `Emf` / `Eo` | 平均场 / KS 能量（来自 `system.get().Eo`） |
| `X` | 精确交换 |
| `SX-X` / `CH` | 屏蔽交换修正与 Coulomb-hole（或全频 residual / 积分） |
| `Sig` | `X + SX-X + CH` |
| `Vxc` | 交换关联参考 |
| `Eqp0` | `Eo + Sig - Vxc` |

---

## 5. 内存结果 `E`

不是文件，而是 MATLAB 主返回值：

```matlab
E = qp.launcher(config);
% 或
E = qp_driver('./SAVE');
```

字段：`Eqp`、`Ex`、`Esx_x`、`Ecoh`、`Eqp0`（空）、`fout`（`qp.dat` 路径）。  
详见 [`class_reference.md`](class_reference.md)。

---

## 6. `r-<prefix>.log`

**路径：** `<output_dir>/r-<prefix>.log`  
**写出：** `service/+output`（在 `input_driver` 中由 `display_input_summary` 打开）

**说明：**  
主运行报告 / 日志（Yambo 风格 OF 命名）。

- 默认 `prefix = 'QP'` → `r-QP.log`
- 示例：`r-Si_gamma.log`、`r-Si2_uc.log`
- 已存在时轮转：`r-Si_gamma_01.log`，…

详细程度由 `CONTROL.log_level`（`0`–`2`）控制。  
`CONTROL.log_file` 为旧的独占日志选项；当前驱动使用的实时报告是 `r-<prefix>.log`。

---

## 7. ISDF 报告与 checkpoint

### 7.1 文本报告

**目录：** `./isdf_report/`（相对 case 工作目录）  
**文件名**（`filename_map`）：

| 模式 | 作用 |
|------|------|
| `o-ISDF_cond` | 条件数 / conditioning 报告 |
| `o-ISDF_HF_id%d` | 指定 ISDF id 的 HF 校验报告 |
| `o-ISDF_adaptive_id%d` | 自适应 ISDF 进度 / 诊断 |

在开启 ISDF（`ISDF.isisdf`）且走对应校验 / 自适应路径时生成。

### 7.2 自适应 checkpoint

**路径：** `<storage_dir>/isdf_adaptive_checkpoint_<desc>_idN.mat`  
（如 `vc` / `vn` / `nn`）

用于恢复昂贵的自适应 ISDF 构造。

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-11 | ZZ | 中文用户版；拆分为 `doc_user_ZH` / `doc_user_EN` |
| 2026-08-08 | Zhengbang | Rewrite for SAVE / qp.dat / r-*.log layout |
