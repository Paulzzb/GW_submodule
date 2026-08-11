# 用户口文档边界（草稿，可改）

日期：2026-08-11  
用途：先定「发布/用户入口能看到哪些 docs」，再改内容与链接。本文件本身可随讨论逐条修改。

## 基本原则（已定）

1. **用户看不到 `doc_devel/`** — 发布包物理排除；用户文档内禁止链接或「见 devel」。
2. 用户口只靠：根 `README` + `examples/` + `doc_user/`。

---

## 1. 用户可见面（已定）

| 路径 | 角色 | 用户是否看到 | 说明 |
|------|------|:------------:|------|
| 根 `README.md` | 项目名片（安装/引用/许可证） | 是 | 入口；只链 `doc_user/` / `examples/` |
| `examples/README.md` | 小展示怎么跑 | 是 | 跟 examples 走 |
| `doc_user/README.md` | 上手：基态、namelist、跑通、看结果 | 是 | **主文档入口** |
| `doc_user/output_description.md` | 产出文件说明 | 是 | |
| `doc_user/GW_input_description.md` | **用户向、简化** namelist | 是 | **以 devel 手册为标准裁剪新建**（不整文件搬家） |
| `doc_devel/GW_input_description.md` | 开发者向完整/复杂参数表 | **否** | 标准全文；发布包不含 |
| `doc_user/class_reference.md` | 极少对外暴露的结构 | 是（极瘦） | **只留 `E` / `qp.dat`**；**不写 `data`**（当前用户侧实质只需 QE 路径） |
| `doc_devel/` 其余 | 内部笔记、架构、完整 class 等 | **否** | 发布包物理排除 |

### 明确不作为用户口

- 整个 `doc_devel/`（发布包排除）
- `tests/`、`test_profile/`、内部 `example_pipe.md`
- 任何「请参阅 doc_devel」式指引（用户树里不允许出现）

---

## 2. 入口关系（示意）

```text
根 README
  ├─→ doc_user/README          （主）
  │     ├─→ output_description
  │     ├─→ GW_input_description   （简化版，仅用户树）
  │     └─→ class_reference        （极瘦）
  └─→ examples/README

doc_devel/                     （发布包无；用户文档零引用）
```

---

## 3. 用户口三份核心叙事

改 `doc_user` 内容时按此收，不先做全局美化：

1. **怎么跑**：`QPstartup` → 准备 QE/KSSOLV → 写 `test` → `input_driver` + `qp.launcher`（及 `examples/run_all`）
2. **怎么填**：`doc_user/GW_input_description.md`（简化 namelist）
3. **怎么看**：`output_description` + 极瘦 `class_reference`（仅 `E` / `qp.dat`；可与 output 合并，不写 `data`）

保留的主跑法示例（与 `examples/run_cohsex` 对齐）：

```matlab
input_driver('./test');
load('./SAVE/config.mat', 'config');
E = qp.launcher(config);
```

两阶段：

```matlab
input_driver('./test');
E = qp_driver('./SAVE');
```

---

## 4. 拍板记录

| # | 议题 | 决定 |
|---|------|------|
| 1 | namelist | **两套**：devel `GW_input_description` 为标准全文；用户侧**按该文裁剪新建** `doc_user/GW_input_description.md`（不整迁）。 |
| 2 | `class_reference` | 用户白名单：**仅 `E`、`qp.dat`**；**不含 `data`**；其余归 devel。 |
| 3 | README §4「扩展框架」 | **不能**写「见 doc_devel」。整节挪 `doc_devel`；用户 README **删除**该节。 |
| 4 | 发布包 | **`doc_devel/` 物理排除**（不只是不链）。 |

§4 微调项已全部拍板，无未决勾选。

---

## 5. 执行状态

1. [x] 用户 namelist 裁剪新建；`class_reference` 仅 E/qp.dat；用户 README 去掉扩展节（→ `doc_devel/EXTENDING.md`）
2. [x] 用户叙事 / 跑法 / 基态说明已改
3. [x] 断链与根 README / examples 入口；**发布包须物理排除 `doc_devel/`**（打包脚本若有，过滤该目录；当前无统一 pack 脚本时作发布约定）
4. [x] 本规划稿已迁到 `doc_devel/USER_DOC_SURFACE.md`（不出现在用户树）

---

## 6. 变更记录

| 日期 | 谁 | 改了什么 |
|------|----|----------|
| 2026-08-11 | AI | 初稿 |
| 2026-08-11 | ZZ+AI | 定：双 namelist；class 极瘦；用户永不看 devel；发布包排除 devel |
| 2026-08-11 | ZZ | namelist：devel 为标准、用户裁剪新建；class 用户只要 E/qp.dat、不要 data |
