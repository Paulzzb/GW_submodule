# 结果结构（用户版）

用户通常只需关心 `qp.launcher` / `qp_driver` 返回的 `E`，以及写出的文本表 `qp.dat`。  
列含义与文件布局详见 [`output_description.md`](output_description.md)。

---

## `E`（`qp.launcher` 返回值）

| 字段 | 说明 |
|------|------|
| `Eqp` | 准粒子能量：`Ex + Esx_x + Ecoh`（**Ry**，选中能带 × k） |
| `Ex` | 精确交换 Σ_x（Ry） |
| `Esx_x` | 屏蔽交换修正（COHSEX 的 SX−X，或全频 residual）（Ry） |
| `Ecoh` | Coulomb-hole / 积分贡献（Ry） |
| `Eqp0` | launcher 不填；写表时用 `Eo + Sig − Vxc` |
| `fout` | 写出的 `qp.dat` 路径（失败时为空） |

```matlab
E = qp.launcher(config);
% 或
E = qp_driver('./SAVE');
```

内部数组为 **Ry**；`qp.dat` 中为 **eV**（`ry2ev ≈ 13.6057`）。

---

## `qp.dat`

- 路径：一般在 `CONTROL.output_dir` 下，文件名固定为 `qp.dat`
- 由 `qp.fout` 写出（`qp.launcher` 内会调用）
- 列包括能带指标与 Emf / Eo / X / SX−X / CH / Sig / Vxc / Eqp0 等

完整表头与单位见 [`output_description.md`](output_description.md#4-qpdat)。

---

### Changing log

| Date | Name | Changes |
|------|------|---------|
| 2026-08-11 | ZZ | 用户版：仅保留 `E` / `qp.dat` |
