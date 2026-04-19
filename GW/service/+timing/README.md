# `+timing` 服务包说明

面向 Yambo `services/timing` 的 MATLAB 移植骨架：持久状态挂在 `timing.manager` 中，类型为 `timing.base.timing_m`。风格尽量兼容 **MATLAB R2008a**（`char`、`addParamValue` 等）。

---

## 顶层入口（与 `symmetry` / `FFT` 等包一致）

| 文件 | 作用 |
|------|------|
| **`driver.m`** | 分配全局/内部 clock 列表（对应 Yambo `timing_allocate`），写入 `manager`。不读 GW `config`；可选 `timing.driver('nclock_max', N)`。 |
| **`get.m`** | `timing.manager('get')` 的简写，取出当前的 `timing_m`。 |
| **`save2mod.m`** | 将 `timing.base.timing_m` 写回 `manager`。 |
| **`free.m`** | 清空 `manager` 中的持久数据。 |
| **`manager.m`** | 内部 `persistent` 容器；支持 `'get'` / `'save2mod'` / `'free'`。 |
| **`report.m`** | 输出时间报告：**三道横线**、中间两段——上段为报告标题 + 列名，下段为各 clock 行与 `TOTAL (listed)`。若参数或 `tm.live.report_logfile` 给出日志路径则 **追加** 写入该文件，否则写入命令行。可选 `'tm', timing_m` 只对给定对象排版、不写回 `manager`。 |

---

## 进度条 / LIVE（Yambo `LIVE_timing*`）

| 文件 | 作用 |
|------|------|
| **`LIVE_timing.m`** | 对外主入口，对应 `LIVE_timing.F`：`()` 关闭；`(steps)` 累加步数；`('标签', totalSteps)` 启动；可选第三参 `depth` 设置 `memory_steps`。 |
| **`timing_get_time.m`** | 对应 `TIMING_get_time.F` 的单进程版：维护 `live.cput_seg` / `cput_sec` / `cput_tot`（均为 **1×2**，列 1 为耗时、列 2 为 `cputime` 锚点）。`varargin` 传入 `'INIT'`、`'INIT_SEG'`、`'SEG'` 等标志。 |
| **`timing_string.m`** | 对应 `TIMING_string.F`：把秒数格式化为 `d/h/m/s` 片段组成的 `char`。 |

### `private/`（仅由包内调用）

| 文件 | 作用 |
|------|------|
| **`live_timing_activate.m`** | 打开 LIVE、重置计数与段计时起点，打第一行进度（`--` 占位）。 |
| **`live_timing_add.m`** | 按步数更新 `steps_done`、hash 档与 ETA；满足 Yambo 规则（hash 前进 + 最小报告间隔等）时刷新一行。 |
| **`live_timing_update.m`** | 拼 `名称 \|####... \| [xxx%] elapsed(E) total(X)` 并 `fprintf`。 |
| **`live_timing_close.m`** | 将步数补满到 `time_steps`，然后关闭 LIVE。 |

---

## `+base/` 数据类型（对应 Yambo `mod_TIMING*`）

| 文件 | 作用 |
|------|------|
| **`clock_m.m`** | 单个计时器槽：`name`、`call_number`、`total_time`、`running` 等（`TYPE(clock)`）。 |
| **`clock_list_m.m`** | 时钟列表：`preallocate` 预分配槽位，`allocate_next_clock` 登记新名字（`TYPE(clock_list)` + `CLOCK_list_allocate` / `CLOCK_allocate`）。 |
| **`timing_m.m`** | 总状态：`global_list`、`internal_list`、`TIMING_verb`、`live`、`logo`、`alloc`；静态方法 `create_default` 完成默认分配。 |
| **`timing_live_m.m`** | LIVE 模块变量：`nhash`、`time_steps`、`steps_done`、`cput_*`、`timing_name`、`live_report_min_seconds`（控制 LIVE 最小刷新间隔，0 表示 hash 变即打行）、`report_logfile`（非空时 `timing.report` 默认追加到该文件）等。 |
| **`timing_logo_m.m`** | Logo 缓冲行（`TIMING_logo`），`logo_line` 为 cell 列。 |

---

## 与 relay / reset

- **`relay.get_relay_config`** 中已注册 `timing.main`，可与其它 service 一并 `collect` / `restore`。
- **`service_reset_persistent`** 会调用 `timing.free()`。

---

## 示例脚本

- `GW/test_profile/test_timing/demo_timing_base.m`：clock 列表与 relay 快照。
- `GW/test_profile/test_timing/demo_live_timing.m`：`LIVE_timing` 条状进度示例。
