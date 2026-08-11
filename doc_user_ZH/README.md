# 快速上手指南

本指南面向**使用者**：如何准备基态、填写 namelist、跑通计算并查看结果。

English version: [`../doc_user_EN/README.md`](../doc_user_EN/README.md)

| 文档 | 内容 |
|------|------|
| [`GW_input_description.md`](GW_input_description.md) | namelist 参数（用户精简版） |
| [`output_description.md`](output_description.md) | 输出文件说明 |
| [`class_reference.md`](class_reference.md) | 结果结构 `E` / `qp.dat` |
| [`../examples/README.md`](../examples/README.md) | 小展示算例（ISDF vs dense） |

---

## 1. 你会用到的目录

```text
QP_root/
├── QPstartup.m         # 启动：加 path、必要时编译 MEX
├── driver/             # input_driver / qp_driver
├── examples/           # 展示算例（推荐先跑这里）
├── doc_user_ZH/        # 本目录：中文用户文档
├── doc_user_EN/        # 英文用户文档
├── interfaces/         # 读 QE / KSSOLV 基态
├── input/              # 解析 namelist
├── service/ · packages/ · util/ · common/ · …
└── …
```

运行时主要产物：

- `SAVE/`（或你设的 `storage_dir`）：`data.mat`、`config.mat`、`relay_stage.mat`
- `qp.dat`、`r-<prefix>.log`
- 可选：`isdf_report/`

---

## 2. 计算流程（简图）

```text
  QE qe.save 或 KSSOLV groundstate.mat
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
         service 构建 → relay_stage.mat
              |
              v
         qp.launcher / qp_driver
              |
              v
         E  +  qp.dat  +  r-*.log
```

`qp.launcher` 按 `FREQUENCY.frequency_dependence` 选择路径：

| 值 | 路径 |
|---:|------|
| `-2` | COHSEX（Gamma）：`gw.x_Gamma` + `gw.cohsex_Gamma` |
| `2` | 全频 CD（Gamma）：`gw.fullfreq_cd_res_Gamma` + `gw.fullfreq_cd_int_Gamma` |

---

## 3. 快速开始

### Step 1：环境初始化

在 MATLAB 中进入仓库根目录：

```matlab
QPstartup;
```

会重置 path，并在需要时编译 ISDF MEX。

### Step 2：准备基态（推荐 QE）

- QE 需支持 **HDF5** 波函数输出。
- 用 `pw.x` 做 `scf`（必要时再 `bands` / `nscf`）。
- 用 `pw2bgw.x`（或等价流程）导出 `vxc.dat`。
- 将下列文件放在同一目录（如 `./qe.save/` 或 `examples/cases/qe.save/`）：

```text
charge-density.hdf5
data-file-schema.xml
wfc*.hdf5
vxc.dat
```

namelist 中设置：

```text
&CONTROL
  groundstate_dir = './qe.save',   % 或 '../qe.save'
  groundstate_type = 'qe',
  storage_dir = './SAVE',
  output_dir = './',
  prefix = 'Si_gamma',
END &CONTROL
```

**KSSOLV**：（暂不支持）目录中需有 `groundstate.mat`（变量名 `groundstate`），并设 `groundstate_type = 'kssolv'`。

参数细节见 [`GW_input_description.md`](GW_input_description.md)。完整可跑示例见 `examples/cases/*/test`。

### Step 3：运行

推荐（与 `examples/run_cohsex.m` 一致）：

```matlab
input_driver('./test');                 % 生成 SAVE/ 并构建 service
load('./SAVE/config.mat', 'config');
E = qp.launcher(config);                % 执行 QP
```

两阶段（可在另一次 MATLAB 会话中从 stage 恢复）：

```matlab
input_driver('./test');
E = qp_driver('./SAVE');
```

只想先看展示对比：

```matlab
cd examples
run_all
```

### Step 4：检查结果

| 文件 | 含义 |
|------|------|
| `r-<prefix>.log` | 运行报告 |
| `qp.dat` | 准粒子能带表（eV） |
| `SAVE/data.mat` | 基态数据 |
| `SAVE/config.mat` | 解析后的配置 |
| `SAVE/relay_stage.mat` | service 缓存 |
| `isdf_report/` | ISDF 诊断（若开启） |

详见 [`output_description.md`](output_description.md)；`E` 字段见 [`class_reference.md`](class_reference.md)。

### Stage 缓存（简要）

- 若已有 `relay_stage.mat` 与 `data.mat`，`input_driver` 会复用 `data`，但 **始终** 按 namelist 重建 `config`。
- 改 CUTOFFS / ISDF 等通常不必删 stage；改基态或对称性相关设置时，请删除 `data.mat` 与 `relay_stage.mat` 后重跑。

---

## 4. 相关文档

| 文档 | 内容 |
|------|------|
| [`GW_input_description.md`](GW_input_description.md) | namelist（用户精简版） |
| [`output_description.md`](output_description.md) | 输出文件 |
| [`class_reference.md`](class_reference.md) | `E` / `qp.dat` |
| [`../examples/README.md`](../examples/README.md) | 展示算例 |

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-11 | ZZ | 拆分为 `doc_user_ZH` / `doc_user_EN` |
| 2026-08-11 | ZZ | 用户口：去掉开发扩展节；链到用户版 namelist / 精简 class |
| 2026-08-08 | Zhengbang | Rewrite for service/packages architecture |
