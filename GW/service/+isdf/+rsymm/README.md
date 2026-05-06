# `+isdf/+rsymm` — 空间对称与 bundle

## `bundle_refresh.m`

在已有 `bundle_struct` 上增加 `N_new` 个采样种子（由 `indices_new` 指定细 FFT 格点），扩展闭包内的 bundle 点与 `R_rot_in_bundle`，并写回 `isdf_data.bundle_struct`。

### 接口

| 形参 | 含义 |
|------|------|
| `id` | 源 ISDF 槽位；从此槽读当前 `bundle_struct` 与采样计数。 |
| `N_new` | 本次新增的 **采样行** 个数（与 `indices_new` 长度一致）。 |
| `indices_new` | 长度 `N_new`；每个元素是 **细 FFT 目录行** 下标，用于 `fft_data.Rgrid_RLU(indices_new(inew), :)` 取种子点的 RLU（与 `nr = prod(fftgrid)` 的格点编号一致，不是 bundle 行号）。 |
| `idnew`（可选） | 缺省或空则 `idnew = id`，原地更新。若 `idnew ~= id`，先把 `id` 的整份 `isdf_m` 拷到 `idnew`（槽空时），再只对 `idnew` 更新 `bundle_struct`。 |

### 尺度与计数

| 变量 | 含义 |
|------|------|
| `nr` | 细 FFT 实空间格点数，`prod(fftgrid)`。 |
| `nsym` | 空间对称操作个数（`symmetry.nsym`）。 |
| `N_sampling_old` / `N_sampling_new` | 刷新前后 **ISDF 采样行总数**（粗 + 细）；`N_sampling_new = N_sampling_old + N_new`。 |
| `N_coarse` | 粗格采样行数；来自 `bundle_struct.N_coarse`（细段为 `N_coarse+1 : N_sampling_new`）。 |
| `Nb_old` / `Nb_new` | 刷新前后 **bundle 上互异格点个数**（`R_grid_bundle` 的行数）。 |

### 主数组（工作区 → 写入 `bundle_struct`）

| 变量 | 形状 / 类型 | 含义 |
|------|-------------|------|
| `Rgrid_b_new` | `nr × 3`（实际只用前 `Nb_new` 行） | 每个 bundle 点的 RLU 坐标（`fft_data.Rgrid_RLU` 中的行，在环面上取整/取模后与 FFT 一致）。 |
| `WF_b_new` | 与 `wf_data.c` 同型，再裁到 `Nb_new` 行 | 各 bundle 点上的 `ψ(ib, ik, ispin)`。 |
| `s2b_new` | `N_sampling_new × 1`，`int32` | 写入 **`bundle_struct.sampling2bundle`**：**采样行 `ialpha` → bundle 行 `irb`**（`1…Nb_new`），**不是**细网格线指标 `1…nr`。 |
| `Rrot_b_new` | 先 `nr × nsym`，保存时裁为 `Nb_new × nsym` | 写入 **`bundle_struct.R_rot_in_bundle`**：`Rrot_b_new(irb, isym)` = 在对称 `isym` 下，bundle 点 `irb` 的像所对应的 **bundle 行号**。 |
| `finegrid2newb` | `nr × 1` | **细网格线指标** `1…nr` → **仅对本轮新加入的 bundle 行**（`Nb_old+1:Nb_new`）建立的逆映射：用于把 `indR`（旋转后的 RLU 对应的线指标）映回 bundle 行号，填满 `Rrot_b_new(Nb_old+1:Nb_new, isym)`。 |
| `irbninfft_i` | 长度 `Nb_new - Nb_old` | 新 bundle 行 RLU 对应的 **细 FFT 线指标**（与 `finegrid2newb` 的键一致）。 |

### 整数 RLU 子集与 `bundle2finegrid`

| 变量 | 含义 |
|------|------|
| `ind_b_in_finegrid` | 长度 `Nb_new` 的缓冲；前 **`Nb_in_finegrid`** 个有效元素为 **bundle 行号 `irb`**，满足 `Rgrid_b_new(irb,:)` 与整数 RLU 一致（`|R - round(R)|` 很小）。**注意：存的是 bundle 下标，不是 `1…nr`。** |
| `Nb_in_finegrid` | 上述 **bundle 行** 的个数；后续旋转校验、与 `WF_apply_symm` 对拍时，只在这些行上做细格提取。 |
| `bundle2finegrid` | 长度 `Nb_in_finegrid`，`int32`；**细 FFT 线指标 `1…nr`**：由 `R_grid_bundle(ind_b_in_finegrid(k), :)` 按 `fftgrid` 公式换算，用于 `wf_dir(bundle2finegrid)` 与 bundle 上 WF 对比。 |

### 轨道扩展用到的临时量

| 变量 | 含义 |
|------|------|
| `r_new` | 当前新种子的 RLU（`single`，来自 `Rgrid_RLU`）。 |
| `tmp` | 若种子已在当前 `Rgrid_b_new(1:Nb_new,:)` 中，为其 **bundle 行号**；否则走轨道分支。 |
| `r_rot` | `nsym×3`；各对称下旋转后的 RLU（再 `mod` 到环面）。 |
| `ind_rot_finegrid` / `ind_orbit_finegrid` | 旋转后线指标及 **unique** 后的轨道点集，用于追加 `Rgrid_b_new` / `WF_b_new` 行。 |
| `new_generator` | 记录新种子的 RLU（调试用 / 扩展用）。 |

### 其它槽位

| 变量 | 含义 |
|------|------|
| `isdf_data1` | `isdf.get(id)`，在 `idnew ~= id` 时仍保留 **源槽 `id`** 的 `coeff_seper`、`R_rot_coarse` 等，供后面粗格对拍。 |

### 与 `gen_bundle.m` 的关系

`gen_bundle` 从零构建首个 `bundle_struct`；`bundle_refresh` 在其上 **增量** 增加采样与闭包点，并保持 **`sampling2bundle` / `R_rot_in_bundle` / `R_grid_bundle` / `WF_bundle`** 的约定一致。
