# QP框架开发者快速上手指南

本指南旨在帮助新开发者快速理解并扩展 QP（准粒子）计算框架。系统采用模块化设计，结构清晰，支持多人协作开发。

---

## 1. 文件夹结构简介

```
QP_root/
├── common/             # 通用功能函数
├── driver_profile/     # 驱动主程序模块
├── input/              # 输入文件解析与默认参数
├── src_profile/        # 核心类与数值方法
├── GW_profile/         # GW计算辅助子模块
├── util_profile/       # 输出与日志函数
├── test_profile/       # 测试脚本
├── example_profile/    # 使用例子
├── QPstartup.m         # 启动脚本（设置路径）
```

---

## 2. 快速开始


### Step 1：环境初始化

在 MATLAB 中进入 `GW` 文件夹内，并执行

```matlab
QPstartup;
```

### Step 2：准备输入文件

#### A. KSSOLV

- **请务必使用根目录下的 `./kssolv` 进行 scf 计算!!!!!!**

A.1 使用 `KSSOLV` (请使用根目录下 `./kssolv` 地址中的版本) 进行基态计算
```matlab
% set mol
% set options_scf
[mol,H,X0,info] = scf(mol, options_scf);
```

A.2 使用根目录中函数`save_groundstate_to_GWformat`保存数据至目录`my_dir`：

```matlab
save_groundstate_to_GWformat(mol, H, X0, info, 'my_dir');
```

A.3 在当前文件夹下，编写 `input` 文件

```text
&CONTROL
  groundstate_dir = 'my_dir',
  groundstate_type = 'kssolv',
END &CONTROL
```

#### B. Quantum ESPRESSO

- 用户使用的 **Quantum ESPRESSO (QE)** 版本应支持 **HDF5 输出**。

2.1. 使用 `pw.x` 进行两步计算：
   - 第一步：`calculation = 'scf'`；
   - 第二步：`calculation = 'bands'`。
   
   其中，输入文件中 **k 点** 的设置应如下所示：

```text
K_POINTS crystal
1
0.0 0.0 0.0 1.0
``````

2.2 使用 `pw2bgw.x` 导出所需的 `vxc.dat` 文件。

2.3 将以下文件（其中前三个由 `pw.x` 生成，默认位于 `$PREFIX.SAVE` 文件夹内，最后一个由 `pw2bgw.x` 生成）
```text
charge-density.hdf5  data-file-schema.xml  wfc1.hdf5  vxc.dat 
```
复制或链接至目标目录 (例如 `./qeinput/`)。


2.4 在当前目录编写 `input` 文件如下：

```text
&CONTROL
  groundstate_dir = './qeinput',
  groundstate_type = 'qe',
END &CONTROL
```


### Step 3：运行主程序

```matlab
input_driver('./input');  % 读取输入配置
qp_driver('./');          % 执行QP计算
```

### Step 4：检查结果

用户可以在 `QPlog` 中查看计算日志，并在 `QPenergy` 中查看计算的能带结果。

---


## 3. 扩展框架功能（开发新模块）

### Step 1：添加控制参数

在以下文件中注册控制参数 `enable_your_module`：

- `input/allowed_parameter_list.m`：

  ```matlab
  def.CONTROL = { ...
    'enable_your_module', ...
  };
  ```

- `input/default_parameter_values.m`：

  ```matlab
  def.CONTROL = struct(...
    'enable_your_module', 0, ...
  );
  ```

- `doc_profile/GW_input_description.md`：添加说明文档（仿照其他参数）

### Step 2：创建模块文件夹并添加路径

- 在 `QP_root` 下新建 `your_module/` 文件夹；
- 将 `.m` 函数放入；
- 修改 `QPstartup.m`，在预留位置加入路径：
  ```matlab
  % add_mpaths_only([CPATH 'your_module/']);
  ```

### Step 3：在 `qp_driver.m` 中挂接模块调用

在预设 Hook 区域添加以下代码：

```matlab
if config.CONTROL.enable_your_module
  QPlog('Your module is enabled. Starting execution...', 1);
  GWenergy = your_kernel(GWinfo, config);
end
```

### Step 4：编写你的模块函数 `your_kernel.m`

模板示例见附录，确保使用统一输入输出规范，并在函数头部调用 `QPlog_push` 管理日志。

---

## 4. 测试与调试建议

- 在 `test_profile/` 下为每个功能创建子文件夹与测试脚本；
- 所有测试入口应在 `testall.m` 中注册并验证；
- 所有函数必须模块化、稳定，不得引入副作用或全局变量。

---

## 5. 协作规范

- 函数命名清晰，不使用缩写；
- 参数输入输出通过结构体传递；
- 保持注释清晰，使用中英文解释逻辑；
- **避免对驱动或核心类直接修改**，如需修改请先与维护者沟通。

---

## 6. 联系我们

如需技术支持、功能申请或希望加入开发，请联系项目维护者。

---

# 附录：模块开发模板说明（your\_kernel.m）

本文件为新开发者提供一个规范模板，用于在现有 QP 框架中开发新模块。

---

## 一、开发流程总览

### 1. 增加控制参数

#### （1）在 `allowed_parameter_list.m` 中添加参数名

位于 `QP_root/input/allowed_parameter_list.m`，在 `CONTROL` 模块下添加：

```matlab
def.CONTROL = { ...
  'enable_your_module', ...  % 其他参数名
};
```

#### （2）在 `default_parameter_values.m` 中添加默认值

位于 `QP_root/input/default_parameter_values.m`，在 `CONTROL` 模块下添加：

```matlab
def.CONTROL = struct(...
  'enable_your_module', 0, ...  % 默认关闭
  % 其他控制参数
);
```

#### （3）在输入说明文档中添加说明

打开 `QP_root/doc_profile/GW_input_description.md`，按 `CONTROL` 模块格式添加文字说明。

> 若需在现有模块中添加控制参数，请参考以上三个步骤。

---

### 2. 增加额外参数模块（可选）

当现有模块无法满足参数需求时，可新增参数块模块 `MYBLOCK`，流程如下：

#### （1）添加参数名

```matlab
def.MYBLOCK = { ...
  'p1', ...  % 示例参数名
};
```

#### （2）添加默认值

```matlab
def.MYBLOCK = struct(...
  'p1', 1  % 默认值
);
```

#### （3）添加说明文档

在 `GW_input_description.md` 中按 `CONTROL` 模块格式添加 `MYBLOCK` 说明。

> 约定规范：模块名称全大写，参数名称全小写。

---

### 3. 创建模块文件夹并加入搜索路径

在根目录 `QP_root/` 下新建模块文件夹：

```
QP_root/mymodule/
```

将你的模块函数放置其中，并在 `QP_root/QPstartup.m` 中添加路径：

```matlab
% add_mpaths_only([CPATH 'mymodule/']);
```

系统已在对应位置预留 Hook，开发者只需取消注释并修改路径即可。

---

### 4. 在 `qp_driver.m` 中插入调用逻辑

我们在 `qp_driver.m` 中为开发者预留了统一的调用接口区：

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Hook] Insert your custom module calls below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

请将你的模块调用逻辑插入其中，建议格式如下：

```matlab
if config.CONTROL.enable_your_module
  QPlog('Your module is enabled. Starting execution...', 1);
  GWenergy = your_kernel(GWinfo, config);
end
```

#### 调用规范：

- 函数应使用 `GWinfo` 和 `config` 为输入，输出为 `GWenergy`
- 若需额外信息，请联系维护者扩展结构体字段

---

## 二、`your_kernel.m` 模板示例

新建模块函数应遵循如下模板结构：

```matlab
function GWenergy = your_kernel(GWinfo, config)
  % your_kernel - 示例模块函数

  cleanup = QPlog_push('your_kernel');  % 入栈日志模块名

  QPlog('Your kernel calculation started...', 1);

  if (config.MYBLOCK.p1 > 0)
    QPlog('p1 is positive, Performing main computation...', 1);

    % 模块主要计算部分
    result = rand(10);

    QPlog('Computation completed.', 2);
  else
    QPlog('p1 <= 0, main calculation disabled.', 1);
    result = [];
    GWerror('p1 <= 0, main calculation disabled.');
  end
end
```

#### 日志建议层级：

- `0`：强制输出（如异常、终点信息）
- `1`：主流程信息（如启动、结束）
- `2`：详细调试信息（中间步骤）

> 注意：必须调用 `QPlog_push` 和 `onCleanup(QPlog_pop)` 进行日志模块管理。

---

## 三、模块结构建议（可选拆分）

若模块结构复杂，建议拆分为多个子函数：

```
your_kernel.m         % 主接口
your_kernel_core.m    % 主要计算逻辑
your_kernel_post.m    % 后处理逻辑
```

建议将这些模块放置于：

```
QP_root/your_module/
```

---

## 四、调用接口示例

你可以在 `qp_driver.m` 中调用如下逻辑来连接模块：

```matlab
if isfield(config, 'YOURMODULE') && config.YOURMODULE.enabled
  QPlog('Your module is enabled. Starting execution...', 1);
  result = your_kernel(GWinfo, config);
end
```

如有进一步需求或疑问，请联系主维护者以统一扩展数据结构与接口。

