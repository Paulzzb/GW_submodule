# 干净回归测试树（方案 S1）

本目录是新的基准测试入口，与杂乱的 `test_profile/` 分离。  
当前只搭 **MATLAB 侧输入/输出结构**；QE `*.save` 由你随后拷入。

## 布局

```text
example/
  README.md
  run_all.m              % 编排 Si_gamma / Si_k（缺 groundstate 时跳过）
  run_cohsex.m           % 共享：QPstartup → input_driver → qp.launcher
  run_sc.m               % SC 演示：Si2_uc → Si8_sc
  cases/
    Si_gamma/            % Si，无 k 点（enable_k_points = false）
    Si_k/                % Si，带 k 点
    Si2_uc/              % SC 原胞源（小 Si 原胞 + adaptive ISDF）
    Si8_sc/              % SC 超胞（formal 2×2×1 + SC_ISDF）
    SrTiO3_k/            % 搁置：QE 输出过大，不随仓库提交/不进 run_all
```

当前默认回归只有 **Si_gamma** 与 **Si_k**。  
SC 演示另跑：`run_sc`（摘自 `test_profile/test_SC`，体积刻意压小）。

每个 case：

```text
<case>/
  README.md
  test                   % namelist（唯一输入配置）
  qe.save/               % 占位：请拷入 QE 输出（见该目录 README）
  SAVE/                  % 运行产物（gitignore，勿提交）
  ref/                   % 可选：日后放参考 qp_*.dat
```

## 约定（S1）

- case 目录 **不放** `test_input.m` / 重复 `QPstartup` 样板
- 启动逻辑只在 `tests/run_*.m`
- 物理参数以各 case 的 `test` 为准；数值（能带数、截断等）在 QE 数据到位后再定

## 用法（数据到位后）

```matlab
cd <repo_root>/example
run_all
% 单 case:
run_cohsex(fullfile(pwd,'cases','Si_k'))
% SC 演示（先原胞再超胞）:
run_sc
```

需先保证仓库根在 MATLAB 路径上，或由 `run_cohsex` / `run_sc` 自行 `QPstartup`。

## 与旧树关系

| 目录 | 角色 |
|---|---|
| `tests/` | **新**回归闸门（本阶段建设） |
| `test_profile/` | 旧实验场；Phase A 按文件夹审查，不作为新 `testall` 主入口 |
| 根目录 `testall.m` | 待 `tests/` 可用后改为调用 `tests/run_all` |
