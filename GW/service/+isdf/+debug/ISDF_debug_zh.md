# ISDF 内联调试控制（`+isdf/+debug`）

（英文版：`ISDF_debug.md`）

## 会话快照（推荐用法）

调试相关开关从 **`config.ISDF` 只读取一次**，并缓存在当前 MATLAB 会话中：

- **`isdf.debug.init_from_config(config)`** 在 **`service_driver`** 启动时调用；在 **`isdf.adaptive.run_adaptiveisdf`** 中于 relay `restore` 之后也会调用（不经过 `service_driver` 的路径）。
- 初始化之后，各调用点使用 **`isdf.debug.on(tag)`** 与 **`isdf.debug.react(...)`**，**不再传入 `config`**。
- **`isdf.debug.clear()`** 清空缓存；已在 **`service_reset_persistent`** 中注册，与其它 ISDF 持久态一并清理。

| 字段（`config.ISDF`） | 默认值 | 含义 |
|----------------------|--------|------|
| `debug_checks` | `false` | 带 tag 的内联检查总开关。 |
| `debug_level` | `'error'` | `'error'` → 违反时 `error()`；`'warn'` / `'warning'` → `warning()`。 |
| `debug_tags` | `[]` | 空：在 `debug_checks` 为真时**所有** tag 生效；非空：仅列表中的 tag（不区分大小写）。 |

默认值在 `GW/input/default_param_values.m` 中设置。

---

## API 说明

### `isdf.debug.init_from_config(config)`

从 `config.ISDF` 中快照 `debug_*` 字段。若无 `ISDF` 字段则安全（检查保持关闭）。

### `isdf.debug.clear()`

清除快照（在下次 `init_from_config` 之前检查均视为关闭）。

### `tf = isdf.debug.on(tag)`

放在**昂贵计算之外**（`if isdf.debug.on(tag) ... end`）。**仅**使用会话缓存。若从未调用 `init_from_config`，恒为 `false`。

### `isdf.debug.react(violation, msg [, tagForId])`

当 `violation` 为真时，按缓存的 `debug_level` 发出 `warning` 或 `error`。可选 `tagForId` 用作告警标识尾部（会规整为 `[a-zA-Z0-9_]`）。

### 重载（可选：不经过 `service_driver` 的脚本）

- **`isdf.debug.on(config, tag)`** — 直接读 `config.ISDF`；**不会**写回会话缓存。
- **`isdf.debug.level(config)`** / **`isdf.debug.level()`** — 从传入的 config 或缓存取等级。
- **`isdf.debug.react(config, violation, msg [, tagForId])`**
- **`isdf.debug.check(tag, fh, msg [, id])`** 或 **`isdf.debug.check(config, tag, fh, msg [, id])`**

### 内部接口

- **`isdf.debug.cache`** — 持久存储；一般业务代码勿直接依赖。

---

## 用法（模式 B + `react`）

```matlab
if isdf.debug.on('rsymm/bundle_refresh')
  % 仅在开启调试时做昂贵检查
  isdf.debug.react(norm(a - b) > tol, 'message', 'my_check_id');
end
```

### 当前使用的 tag

| Tag | 位置 |
|-----|------|
| `coeff/coeff_coarse_wf_extract` | `coeff_coarse_wf_extract.m` — 重合粗/细格波函数一致性 |
| `rsymm/bundle_refresh` | `bundle_refresh.m` — 同一开关内：环面旋转、`R_rot` 离散映射、全带 WF bundle/采样、`ib=3:5` 采样与旋转路径对比，以及 **`save2mod` 之后** 对 `bs_new` 的检验（`ib=3:5`、`S_q` 细/粗格与直接法对比）（关闭调试时整段不执行） |

---

## 未走 `service_driver` 的脚本

在拿到与传给 `service_driver` **相同**的 `config` 结构体之后，自行调用一次 **`isdf.debug.init_from_config(config)`**（例如 `load(..., 'config')` 之后）。

---

**英文版路径：** `GW/service/+isdf/+debug/ISDF_debug.md`  
**本文路径：** `GW/service/+isdf/+debug/ISDF_debug_zh.md`  
**更新日期：** 2026-05-02
