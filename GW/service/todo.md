# Service 状态统一 TODO

## 目标

为 service 模块建立一个统一的运行时状态入口和持久化入口，
并在迁移过程中不破坏现有各模块 manager 的行为。

目标流程：

1. 全部 driver 跑一遍。
2. 一条命令保存当前状态。
3. 一条命令重新加载当前状态。

## 命名决策

- 首选包名：`+state`
- 可选别名（强调持久化语义）：`+dbio`

原因：

- `state` 同时覆盖运行时缓存与存储语义，含义更完整。
- `dbio` 更像“只做落盘”，但我们还需要运行时读写。

## 范围（第一阶段）

第一批支持的状态 key：

- `symmetry.main`
- `FFT.main`
- `lattice.r_lat`
- `lattice.d_lat`
- `lattice.k`
- `lattice.q`
- `coulomb.main`（占位，待数据路径就绪后接入）

## 计划

- [ ] Step 0：冻结统一仓库 API 契约。
- [ ] Step 1：创建 `GW/service/+state` 的骨架文件。
- [ ] Step 2：实现内存仓库 API（`get`、`set`、`has`、`clear`、`list`）。
- [ ] Step 3：实现持久化 API（`save`、`load`）并落到 MAT 快照。
- [ ] Step 4：增加桥接函数，从现有 managers 收集数据。
- [ ] Step 5：增加桥接函数，将数据回灌到现有 managers。
- [ ] Step 6：增加一个顶层 helper，执行 `collect -> save`。
- [ ] Step 7：增加一个顶层 helper，执行 `load -> restore`。
- [ ] Step 8：按当前 driver 顺序做 smoke test。
- [ ] Step 9：在 `GW/service/README` 或新文档中补充使用示例。

## API 草案（第一版）

- `state.set(key, value)`
- `state.get(key)`
- `state.get(key, defaultValue)`
- `state.has(key)`
- `state.clear()`
- `state.clear(key)`
- `state.list()`
- `state.save(filePath)`
- `state.load(filePath)`
- `state.collect_from_managers()`
- `state.restore_to_managers()`

## 迁移规则

不要一次性批量替换全部 manager 调用。

按以下顺序推进：

1. 保持当前 manager 逻辑全部可用。
2. 先加桥接层并验证快照流程。
3. 逐模块迁移读取入口。
4. 仅在验证完成后删除重复缓存路径。

## 风险与防护

- 风险：代码变更后类对象不兼容。
  - 防护：快照加入版本字段，load 时做基础 schema 检查。
- 风险：部分加载导致运行态不一致。
  - 防护：load 失败即整体失败，restore 使用全有或全无策略。
- 风险：key 命名漂移。
  - 防护：集中定义 key 常量文件。

## 完成标准

- 全量 driver 运行后，一条命令保存可用。
- 在全新 MATLAB 会话中，一条命令加载并恢复可用。
- 现有模块 manager 行为保持向后兼容。
- symmetry + FFT + lattice 路径的基础 smoke test 通过。
