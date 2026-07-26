# FORMAL benchmark — code changes log

形式复杂度基准：单胞 xml + 超胞缩放 + 随机波函数 / 单调 `Eo`，**不考虑晶体对称性**（`nsym=1`），`frequency_dependence = -2` 走 Gamma 优化路径。

## 新增文件

| 文件 | 作用 |
|------|------|
| `input/load_formal_groundstate.m` | 读单胞 xml，缩放超胞，合成 `psig`/`ev`/`occupation`，强制平凡对称性 |
| `test_profile/test_formal/run_formal_bench.m` | 运行 `input_driver` + `qp_cohsex`，写 `formal_bench_report.txt` |
| `test_profile/test_formal/Si8_formal/test` | 烟雾（Si 2×2×1，32 bands） |
| `test_profile/test_formal/Si16_formal/test` | Si 2×2×2，64 bands（对标 `Si16_from_Si2`） |
| `test_profile/test_formal/LiH_444_formal/test` | LiH 4×4×4，1024 bands（对标 `LiH_444_from_111`） |
| `test_profile/test_formal/LiH_666_formal/test` | LiH 6×6×6 大规模模板（3456 bands） |
| `test_profile/test_formal/FORMAL_BENCH.md` | 本变更记录 |

## 修改文件（最小侵入）

| 文件 | 变更 |
|------|------|
| `input/load_groundstate_info.m` | 增加 `case 'formal'`，可选第三参数 `config` |
| `driver_profile/input_driver.m` | `load_groundstate_info(..., config)` |
| `input/default_param_values.m` | 新增 `&FORMAL` 默认值 |
| `input/allowed_param_list.m` | 登记 `&FORMAL` 参数 |
| `input/set_default_param_value.m` | `formal` 模式校验：`freq=-2`、`use_sc_isdf`、关 `validate_hf`、开 `exact_ch` |
| `input/display_input_summary.m` | 打印 `FORMAL` 块摘要 |
| `service/+isdf/driver.m` | `local_dispatch_gen_tildeVq`：`freq=-2` 时调用 `gen_tildeVq_Gamma` |
| `service/gw_x_Gamma.m` | 函数名修正为 `gw_x_Gamma`（与 `qp_cohsex` 一致） |
| `service/+isdf/cohsex_resolve_ids.m` | 支持 `vc_sc_*` / `vn_sc_*` 等 SC_ISDF 描述符前缀匹配 |
| `service/+isdf/attach_gamma_trunc_factors.m` | `gen_tildeVq_Gamma` 后补全恒等 `CCHq_trunc_factors` |

**未修改**：`+adaptive*`、`SC_ISDF.m`、`gen_tildeVq_Gamma.m` 内部算法；现有 `qe`/`test_SC` 路径行为不变。

## 输入约定

```fortran
&CONTROL
  groundstate_type = 'formal'    ! 单胞 LiH.save（仅 xml）
  enable_k_points = .false.
&FREQUENCY
  frequency_dependence = -2      ! Gamma: gw_x_Gamma + gw_cohsex_test_Gamma
&SUPERCELL
  use_sc_isdf = .true.
  k1/k2/k3 = 扩胞比
  isdf_source_dir = '.../LiH_111/SAVE'
&FORMAL
  wf_seed, eo_e0, eo_delta      ! 随机 wf；Eo(ib)=eo_e0+(ib-1)*eo_delta (Ry)
  run_qp_cohsex = .true.
&ISDF
  validate_hf = .false.
&COHSEX
  exact_ch = .true.
```

## 计算流程

```
load_formal_groundstate (nsym=1, Gamma k)
  → service_driver (FFT, wf, SC_ISDF, …)
  → isdf.driver → gen_tildeVq_Gamma   [if freq=-2]
  → qp_cohsex → gw_x_Gamma + gw_cohsex_test_Gamma
```

## 运行

```bash
cd GW_double_test
matlab -batch "addpath('test_profile/test_formal'); run_formal_bench('test_profile/test_formal/Si16_formal');"
matlab -batch "addpath('test_profile/test_formal'); run_formal_bench('test_profile/test_formal/LiH_444_formal');"
```

`&FORMAL enable_profile = .true.` 或第二参数 `run_formal_bench(case_dir, true)` 会生成 `profile_service/`、`profile_qp_cohsex/`（MATLAB `profsave` HTML），并自动打包 CSS、生成 `index.html` 与 `profile_portal.html`，可在本机浏览器打开。

**本机查看 profiling：**

1. 把整个 case 目录（如 `Si16_formal/`）拷到本机。
2. 用浏览器打开 `profile_portal.html`，或直接进入 `profile_service/index.html` / `profile_qp_cohsex/index.html`。
3. 已有旧 profile 可修复：`matlab -batch "addpath('test_profile/test_formal'); formal_profile_browser_fix('.../profile_service');"`

另会生成 `profile_summary.txt`（按 TotalTime 排序的 Top 函数），可在编辑器里直接读。

**依赖**

- Si 烟雾：`test_Si/Si.save`（xml）、`test_SC/Si2/SAVE`（vc+vn checkpoint）
- LiH：`Materialtest/systems/LiH/1_1_1/LiH.save`；`test_SC/LiH_111/SAVE` 须同时含 `isdf_adaptive_checkpoint_vc_*.mat` 与 `vn_*.mat`（`compute_vc=true` 跑一次 LiH_111）

## 设计说明

- **对称性**：不读取 QE `nsym`；`data.syms` 恒为 `nsym=1`、恒等旋转；`WF_apply_symm` 始终 `isymm=1`。
- **物理量**：结果无物理意义，仅用于 wall-time / 规模验证。
- **与 test_SC 隔离**：`groundstate_type` 分支独立；仅 `isdf.driver` 增加 1 个 `freq=-2` 分发函数。
