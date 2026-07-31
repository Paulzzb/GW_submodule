# SrTiO3_k — 搁置（不进默认回归）

| 项 | 值 |
|---|---|
| 体系 | SrTiO3 |
| k 点 | 是（`enable_k_points = .true.`） |
| 状态 | **deferred**：QE `*.save` 体积大，不适合随仓库提交测试数据 |
| 入口 | 本地自备数据后：`run_cohsex(<此目录>)`（**不**由 `run_all` 调用） |

## 本地使用（可选）

1. 将 QE 输出放到本目录 `qe.save/`（该路径已在 `.gitignore` 中忽略内容）  
2. 按体系修正 `test` 中的能带窗口 / 截断  
3. 单独运行：`run_cohsex(fullfile(pwd,'tests','cases','SrTiO3_k'))`

待有轻量替代（裁剪 save、外部数据仓、或 CI 缓存）后再重新纳入 `run_all`。
