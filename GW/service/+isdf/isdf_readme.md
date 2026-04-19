# ISDF 包 (`+isdf`) — 接口说明

面向 GW 服务的 ISDF 相关工具与粗网格系数路径。实现拆成多个顶层 `.m` 文件，**不使用嵌套子函数**。

---

## 入口与分派

### `isdf.gen_coeff_coarse`

统一粗网格相关操作的入口（文件名 `gen_coeff_coarse.m`）。


| 调用                                                                                  | 返回值                          | 说明                                                                                                                                                                                                             |
| ----------------------------------------------------------------------------------- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `[Nmu, ind_mu] = isdf.gen_coeff_coarse()`                                           | `Nmu` int32，`ind_mu` int32 列 | 在**细网格 FFT 的 R 点列表**上按 gcd/step 规则选子集；阈值 `Nmu_tmp = isdf_ratio * nrep`（`pair_symmetry.nrep`）。实现见 `coeff_coarse_rgrid_indices.m`。                                                                               |
| `[fi, fsz, nd, fc, Nmu, Rc, Rrot] = isdf.gen_coeff_coarse('fft_grid')`              | 7 个输出                        | **规则粗 FFT 子格**：`fftgrid_c` 为整数三维尺寸，`prod(fftgrid_c) > isdf_ratio * nb`（`wave_functions.nb`）；`R_coarse_RLU` 为粗格上 RLU 采样点；`R_rot_coarse` 为对称在粗格上的指标映射。实现见 `coeff_coarse_fft_grid_build.m`、`isdf_coarse_R_rot.m`。 |
| `wf = isdf.gen_coeff_coarse('wf', fftgrid_i, fft_sz, fftgrid_c, Nmu, R_coarse_RLU)` | `wf_on_coarse`               | 将各 IBZ k、带、自旋的波函数映射到粗 RLU 格点（细格实空间 → FFT → `G_table` → 稠密相位矩阵到粗格）。实现见 `coeff_coarse_wf_extract.m`。                                                                                                             |


字符 op 区分大小写不敏感（内部 `lower`）。

### `isdf.gen_indices_coarse`

与 `**isdf.gen_coeff_coarse()`**（无参数）等价，保留给旧脚本（如 `run_coarse_selection_test`）。

---

## 底层实现文件（可直接调用）


| 函数文件                            | 作用                                                                                      |
| ------------------------------- | --------------------------------------------------------------------------------------- |
| `coeff_coarse_rgrid_indices.m`  | R 网格子集 + `ind_mu`；内部用 `divisors_int32`。                                                 |
| `divisors_int32.m`              | 正整数 `n` 的所有正因子（`int32` 列）。                                                              |
| `coeff_coarse_fft_grid_build.m` | 由 `FFT` / `isdf` / `wave_functions` 读服务态，构造粗 `fftgrid_c`、`R_coarse_RLU`、`R_rot_coarse`。 |
| `coeff_coarse_wf_extract.m`     | 粗格波函数抽取（依赖 `lattice`、`FFT`、`wave_functions`）。                                           |
| `isdf_coarse_R_rot.m`           | 给定 `fftgrid_c`，计算 `R_rot_coarse(nr, nsym)`。                                             |
| `isdf_apply_symm_on_coarse.m`   | 在粗格值上施加对称：按 `R_rot_coarse` 置换，时间反演扇区取共轭。                                                |


---

## Manager 多副本池（`N_MAX = 10`，按 **`id`** 索引）

与 Yambo `FFT(FFT_N_max)` 类似：`isdf.manager` 在 **persistent** 里维护 **10 份** `isdf_m`（`id = 1..10`；无并行隔离，仅一个 **当前 `current_id`**）。

| 调用 | 作用 |
|------|------|
| `isdf.get()` / `isdf.save2mod(obj)` | 读/写 **当前 `current_id`** 对应的数据 |
| `isdf.get(id)` | 读指定 `id`，**不**改变 `current_id` |
| `isdf.save2mod(obj, id)` | 写到指定 `id`（省略第二参则写当前 `id`） |
| `id = isdf.isdf_add(desc)` | 找第一个空位，放入占位 `isdf_m`，`desc` 为 `string`，设为当前 `id`，返回该 **`id`** |
| `isdf.isdf_select(id)` | 将 `current_id` 设为 `id`（该 `id` 必须已有数据） |
| `isdf.isdf_current()` | 当前 `current_id`（`int32`） |
| `isdf.isdf_list()` | 每行：池内索引 `id` / `desc` / `assigned` / `empty` |
| `isdf.isdf_nmax()` | 常量 `10` |
| `isdf.free()` | 清空整个池并重置 `current_id` 为 1 |
| `isdf.free(id)` | 只清空 `id`；若删的是当前 `id`，则 `current_id` 改为第一个非空 `id`（否则为 1） |

`isdf.base.isdf_m` 由 manager 写入 **`id`**（池内编号）与 **`desc`**（说明）。`isdf.driver` 若 `config.ISDF` 含 **`desc`** 或 **`note`** 则写入 `desc`；否则默认 `"driver"`，并 `save2mod` 到 **当前 `id`**（通常首次为 `id == 1`）。

`relay.collect` / `restore` 仍针对 **`isdf.get()` 当前 `id`** 的单对象快照。

---

## 与其它模块的关系

- `**isdf.gen_coeff**`（`gen_coeff.m`）中 `method == 'coarse'` 时调用 `**isdf.gen_coeff_coarse()**`（无参），得到 `[Nmu, ind_mu]`。
- `**isdftest/gen_indices_coarse_test.m**` 使用 `**isdf.gen_coeff_coarse('fft_grid')**` 与 `**isdf.gen_coeff_coarse('wf', ...)**`；`isdf_build_tildeVq` / `isdf_coarse_validate_energies` 中粗格对称应调用 `**isdf.isdf_apply_symm_on_coarse**`。
- 依赖服务：`FFT.get`、`isdf.get`、`pair_symmetry.get`、`wave_functions.get`、`symmetry.get`、`lattice.manager`、`system`（若下游需要）等；调用前需已按工程流程初始化（如 `service_driver` / `test_input`）。

---

## 测试

- 粗格流水线演示：`GW/service/isdftest/gen_indices_coarse_test.m`（每次完整重算；中间结果持久化可在数据库层实现）。
- R 网格粗选沙箱：`GW/service/isdftest/run_coarse_selection_test.m` → `isdf.gen_indices_coarse()`。

