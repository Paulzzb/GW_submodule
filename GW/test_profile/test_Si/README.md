# test_Si 测试说明

本目录用于在 Si 数据集上做 `GW`/`ISDF` 相关流程与功能验证。  
下面按“脚本是干什么的”进行简要说明。

## 运行前准备

- 在 `test_Si` 目录下运行。
- 需要已有输入与快照文件：
  - `test`（输入模板）
  - `SAVE/config.mat`、`SAVE/GWinput.mat`
  - `test_relay_stage.mat`
- 若缺少上述文件，先跑一次 `input_driver('./test')` 生成。

## 主脚本说明

- `test_input.m`
  - 读取/检查输入，生成测试所需的 `SAVE` 数据（依赖上层流程）。
- `test_stage.m`
  - 从 `test_relay_stage.mat` 恢复 relay stage 到各 service manager。
  - 主要用于把测试环境恢复到可复现实验状态。
- `run_isdf_driver_smoke.m`
  - `ISDF` 驱动冒烟测试。
  - 执行 stage restore 后调用 `isdf.driver([], config)`，确认流程能跑通。
- `demo_adaptive_weight.m`
  - 演示/检查 adaptive weight 相关流程。
  - 包含 `isdf_exclude_point`、`isdf_schur_update`、`adaptive_weight` 的基本调用链。
- `validate_rgrid_symm_orbits.m`
  - 新增验证脚本：验证两类实空间点集的对称 orbit 构建是否正确。
  - 调用：
    - `rgrid_symm_orbit()`
    - `isdf_r_sampling_symm_orbit(id_nn)`
  - 检查项包括：覆盖完整性、代表元维度、`ir2rep/ir2rot` 索引范围等，并打印简要报告。
- `demo.m`
  - pair symmetry 相关 demo 入口，侧重 `k1-k2` 映射构建与保存（与本次 R-grid orbit 验证不同模块）。
- `test_gw_x_k.m`
  - 调用 `gw_x_k_packages(...)` 做交换项相关测试。

## 子目录说明

- `validation_tests/`
  - 对称性与乘积一致性验证相关脚本：
    - `demo_validation_symmetry_product.m`
    - `val_init_from_stage.m`
  - 以及验证设计文档：
    - `validation_plan.md`
- `Si.save/`
  - QE 导出的 Si 基础数据文件（电荷密度、波函数、赝势等）。
- `SAVE/`
  - 测试中间文件与配置缓存。

## 推荐执行顺序（快速验证）

1. `run_isdf_driver_smoke`
2. `validate_rgrid_symm_orbits`
3. （可选）`demo_adaptive_weight`
4. （按需）`test_gw_x_k` 或 `validation_tests/*`

