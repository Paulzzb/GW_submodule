# test_SC — Supercell ISDF (SC_ISDF) 测试

本目录用于验证 **原胞 adaptive ISDF → 超胞 SC_ISDF 复制 → HF 能量校验** 的完整流程。

- **算法与代码对照文档（LaTeX）：** [`SC_ISDF_report.tex`](SC_ISDF_report.tex)
- QE 结构：**Si** 在 `systems/Si_SCtest/`；**LiH** 在 `Materialtest/systems/LiH/`
- `energy_band_index_max` = **2 × 占据态 band 数**（见下表）

---

## 目录结构

```
test_SC/
├── README.md
├── SC_ISDF_report.tex     # 算法 ↔ 代码 LaTeX 报告
├── run_sc_test.m
├── Si2/                   # Si 原胞 (max band = 8)
├── Si8_from_Si2/          # Si 2×2×1 超胞 (max = 32)
├── Si16_from_Si2/         # Si 2×2×2 超胞 (max = 64)
├── LiH_111/               # LiH 1×1×1 原胞 (max = 16)
├── LiH_222_from_111/      # LiH 2×2×2 (max = 128)
├── LiH_444_from_111/      # LiH 4×4×4 (max = 1024)
└── LiH_666_from_111/      # LiH 6×6×6 (max = 3456，需先完成 QE)
```

---

## Band 范围设置

| Case | $N_{\mathrm{elec}}$ | occupied | `energy_band_index_max` |
|------|----------------------:|---------:|------------------------:|
| Si2 | 8 | 4 | **8** |
| Si8 | 32 | 16 | **32** |
| Si16 | 64 | 32 | **64** |
| LiH_111 | 16 | 8 | **16** |
| LiH_222 | 128 | 64 | **128** |
| LiH_444 | 1024 | 512 | **1024** |
| LiH_666 | 3456 | 1728 | **3456** |

---

## Si 测试（已重跑 2026-06-11）

```bash
cd ~/GW_double_test
matlab -batch "addpath('test_profile/test_SC'); run_sc_test('test_profile/test_SC/Si2');"
matlab -batch "addpath('test_profile/test_SC'); run_sc_test('test_profile/test_SC/Si8_from_Si2');"
matlab -batch "addpath('test_profile/test_SC'); run_sc_test('test_profile/test_SC/Si16_from_Si2');"
```

| Case | nisdf | HF 报告 | 备注 |
|------|-------|---------|------|
| Si2 | 26 | `Si2/isdf_validate_HF_id2.txt` | adaptive，machine precision |
| Si8 | 104 | `Si8_from_Si2/isdf_validate_HF_id1.txt` | SC_ISDF，精度待优化 |
| Si16 | — | `Si16_from_Si2/isdf_validate_HF_id1.txt` | 同上 |

---

## LiH 测试

**依赖：** 先完成 `Materialtest/systems/LiH/1_1_1/` 的 QE（scf → nscf → pw2bgw），`vxc.dat` 在 `LiH.save/` 内。

**Slurm 提交（推荐 largememory 节点）：**

```bash
cd test_profile/test_SC/LiH_111 && sbatch s_test_LiH_111
# 111 完成后：
cd ../LiH_222_from_111 && sbatch s_test_LiH_222
cd ../LiH_444_from_111 && sbatch s_test_LiH_444
# 666 需先跑完 Materialtest/systems/LiH/6_6_6/ 的 QE：
cd ../LiH_666_from_111 && sbatch s_test_LiH_666
```

**groundstate 路径：**

| test 目录 | groundstate |
|-----------|-------------|
| LiH_111 | `Materialtest/systems/LiH/1_1_1/LiH.save` |
| LiH_222 | `Materialtest/systems/LiH/2_2_2/LiH.save` |
| LiH_444 | `Materialtest/systems/LiH/4_4_4/LiH.save` |
| LiH_666 | `Materialtest/systems/LiH/6_6_6/LiH.save` |

**SC 倍率：** 222 → `[2,2,2]`；444 → `[4,4,4]`；666 → `[6,6,6]`；`isdf_source_dir = '../LiH_111/SAVE'`。

> LiH_444 / LiH_666 的 HF 校验 band 数极大（1024 / 3456），运行时间与内存需求显著高于 Si；建议在 largememory 节点分批监控。

---

## SC_ISDF 数据流

```
SC_ISDF → gen_coeff_from_fine_grid → gen_bundle → gen_tildeVq → validate_hf
```

详见 `SC_ISDF_report.tex` 与 [`README.md`](README.md) 历史 Si 章节。

---

## 代码改动索引

| 文件 | 说明 |
|------|------|
| `service/+isdf/SC_ISDF.m` | 新增 |
| `service/+isdf/+coeff/gen_coeff_from_fine_grid.m` | 新增 |
| `service/+isdf/driver.m` | SC 分支 |
| `service/+isdf/+rsymm/gen_bundle.m` | 超胞 scheme |
| `input/*` SUPERCELL 块 | 新增 |
