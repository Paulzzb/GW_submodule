# QP 框架快速上手指南

本指南面向用户与开发者，说明当前版本（service + packages）的目录结构、运行流程与扩展方式。  
参数手册见 [`../doc/GW_input_description.md`](../doc/GW_input_description.md)（该文档已与现版 namelist 对齐）。  
类与数据结构见 [`class_reference.md`](class_reference.md)；输出文件见 [`output_description.md`](output_description.md)。

---

## 1. 文件夹结构简介

```
QP_root/
├── common/             # 通用数值工具（FFT、频率网格、BZ 辅助等）
├── driver/             # 用户入口：input_driver / qp_driver
├── input/              # namelist 解析、默认值、基态加载
├── src/                # 遗留核心（目前主要为 @gvec）
├── service/            # service 层（system / lattice / coulomb / ISDF / output …）
├── packages/           # 物理引擎（+qp, +gw）
├── util/               # filename_map / constant_map
├── example/            # 干净回归用例与 run_*.m
├── test_profile/       # 旧实验场（非新主入口）
├── doc/                # 文档（含 GW_input_description）
├── doc_rev/            # 本目录：按现版代码重写的说明
└── QPstartup.m         # 启动脚本（路径 + ISDF MEX）
```

运行时核心数据：

- 磁盘：`data.mat` + `config.mat` + `relay_stage.mat`
- 内存：service 层 manager（`system.get()`、`lattice.manager(...)` 等）
- 结果：结构体 `E`（由 `qp.launcher` 返回）+ 文本 `qp.dat`

---

## 2. Workflow

```
                 +-------------------------------------------+
                 |  Stage I: Ground-State Data Conversion    |
                 +-------------------------------------------+

    +---------------------------+     +-----------------------------+
    |  Ground-State from        |     |  Ground-State from          |
    |  KSSOLV (暂不支持)         |     |  Quantum ESPRESSO / formal  |
    +---------------------------+     +-----------------------------+
               |                                  |
               x                                  v
    +--------------------------------+    +---------------------------+
    | save_groundstate_to_GWformat   |    | qe.save / formal XML      |
    | → groundstate.mat              |    +---------------------------+
    +--------------------------------+

                 |
                 v

                 +-------------------------------------------+
                 |  Stage II: Input + Service Construction   |
                 +-------------------------------------------+

      +-----------------------------+
      |  namelist (e.g. ./test)     |
      +-----------------------------+
                 |
                 v
      +-----------------------------+
      |     input_driver.m          |
      +-----------------------------+
            /                \
           v                  v
+-------------------+    +--------------------+
|  data.mat         |    |  config.mat        |
|  (groundstate)    |    |  (namelist+defaults)|
+-------------------+    +--------------------+
                 |
                 v
      +-----------------------------+
      |     service_driver.m        |
      +-----------------------------+
                 |
                 v
      +-----------------------------+
      |  relay_stage.mat            |
      |  (+ managers in memory)     |
      +-----------------------------+

                 |
                 v

                 +-------------------------------------------+
                 |  Stage III: Quasiparticle Calculation     |
                 +-------------------------------------------+

      +-----------------------------+
      |  qp_driver / qp.launcher    |
      +-----------------------------+
                 |
                 v
      +-----------------------------------+
      |  gw.x / gw.cohsex_* /             |
      |  gw.fullfreq_cd_*                 |
      +-----------------------------------+
                 |
                 v
      +-----------------------------+
      |  E struct + qp.dat          |
      |  r-<prefix>.log             |
      +-----------------------------+
```

`qp.launcher` 按 `FREQUENCY.frequency_dependence` 路由：

| 值 | 路径 |
|---:|------|
| `-2` | `gw.x_Gamma` + `gw.cohsex_Gamma` |
| `-1` | `gw.x` + `gw.cohsex_multi_k` |
| `2`  | `gw.fullfreq_cd_res_Gamma` + `gw.fullfreq_cd_int_Gamma` |

---

## 3. 快速开始

### Step 1：环境初始化

在 MATLAB 中进入仓库根目录并执行：

```matlab
QPstartup;
```

该脚本会重置 MATLAB path、加入各模块路径，并在需要时编译 ISDF MEX。

### Step 2：准备基态与 namelist

#### A. KSSOLV（暂不支持）

当前发布版不支持以 KSSOLV 作为基态输入；请使用 Quantum ESPRESSO 路径。

#### B. Quantum ESPRESSO

- QE 版本需支持 **HDF5** 波函数输出。

B.1 用 `pw.x` 做 `scf`（必要时再做 `bands`）。

B.2 用 `pw2bgw.x` 导出 `vxc.dat`。

B.3 将下列文件放入目标目录（如 `./qe.save/`）：

```text
charge-density.hdf5
data-file-schema.xml
wfc*.hdf5
vxc.dat
```

B.4 namelist：

```text
&CONTROL
  groundstate_dir = './qe.save',
  groundstate_type = 'qe',
  storage_dir = './SAVE',
  output_dir = './',
  prefix = 'Si_gamma',
END &CONTROL
```

完整参数说明见 `doc/GW_input_description.md`。仓库内可参考 `example/cases/*/test`。

### Step 3：运行主程序

推荐方式（与 `example/run_cohsex.m` 一致）：

```matlab
input_driver('./test');                 % 生成 SAVE/ 并构建 service
load('./SAVE/config.mat', 'config');
E = qp.launcher(config);                % 执行 QP
```

或使用两阶段入口（可在另一次 MATLAB 会话中重启）：

```matlab
input_driver('./test');
E = qp_driver('./SAVE');                % 从 relay_stage 恢复后调用 qp.launcher
```


### Step 4：检查结果

| 文件 | 含义 |
|------|------|
| `r-<prefix>.log` | 运行报告（默认 `r-QP.log`） |
| `qp.dat` | 准粒子能带表（eV） |
| `SAVE/data.mat` | 基态结构体 |
| `SAVE/config.mat` | 解析后的配置 |
| `SAVE/relay_stage.mat` | 缓存的 service 对象 |
| `isdf_report/` | ISDF 诊断报告（若启用） |

详见 [`output_description.md`](output_description.md)。

### Stage 缓存行为（简要）

- 若 `storage_dir` 下同时存在 `relay_stage.mat` 与 `data.mat`，`input_driver` 会复用 `data`，但 **始终** 按 namelist 重建 `config`。
- `service_driver` 对 **lattice / coulomb / ISDF** 每次按当前 config 重建；对 **system / symmetry / FFT / wave_functions** 可从 stage 恢复。
- 修改 CUTOFFS、ISDF 等通常不必删 stage；若基态或对称性相关数据变更，请删除 `data.mat` 与 `relay_stage.mat` 后重跑。

---

## 4. 扩展框架功能（开发新模块）

### Step 1：添加控制参数

在下列文件中注册，例如 `enable_your_module`：

- `input/allowed_param_list.m`
- `input/default_param_values.m`
- `doc/GW_input_description.md`（用户说明）

约定：块名全大写，参数名全小写。

### Step 2：创建模块并加入路径

在根目录新建文件夹，并在 `QPstartup.m` 预留位置加入：

```matlab
% add_mpaths_only([CPATH 'your_module/']);
```

### Step 3：在计算路径中挂接

物理主路径在 `packages/+qp/launcher.m`。若只需旁路模块，可在 `qp_driver.m` 的调用前后，或 `service_driver.m` 中按 config 开关调用。建议接口：

```matlab
if config.CONTROL.enable_your_module
  output.msg('v0s', '%s', 'Your module is enabled.');
  result = your_kernel(config);
end
```

新代码应通过 service manager 取数据（如 `system.get()`、`wave_functions.get()`）。

### Step 4：日志

使用 `service/+output`：

```matlab
cleanup = output.push('your_kernel');
output.msg('v1s', '%s', 'started');
% ...
```

建议层级：`v0` 主流程，`v1` 常规信息，`v2` 调试细节。

---

## 5. 测试与协作

- 新回归入口：`example/`（`run_all` / `run_cohsex` / `run_sc`）
- 旧实验：`test_profile/`（不作为新闸门）
- 函数命名清晰；通过 `config` / service 传递状态，避免隐式全局副作用
- **避免随意改 driver 与核心 service API**；需要变更请先与维护者沟通

---

## 6. 相关文档

| 文档 | 内容 |
|------|------|
| [`../doc/GW_input_description.md`](../doc/GW_input_description.md) | namelist 参数手册（现版） |
| [`class_reference.md`](class_reference.md) | `data` / service 类 / `@gvec` / `E` |
| [`output_description.md`](output_description.md) | `qp.dat`、日志、SAVE 文件 |

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-08 | Zhengbang | Rewrite for service/packages architecture |
