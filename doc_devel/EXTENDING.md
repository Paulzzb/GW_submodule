# 扩展框架（开发者）

从原 `doc_user` 开发扩展节迁出。用户发布包**不包含**本目录。

## 添加控制参数

在下列文件中注册，例如 `enable_your_module`：

- `util/allowed_param_list.m`
- `util/default_param_values.m`
- `doc_devel/GW_input_description.md`（完整参数手册）

约定：块名全大写，参数名全小写。

## 新增基态来源（interface）

在 `interfaces/<source>/` 增加适配器，并在 `interfaces/load_groundstate_info.m` 增加 `case`。  
勿把外部格式解析写进 `service/`；`service` 只消费内部 `data` 契约。

## 创建模块并加入路径

在根目录新建文件夹，并在 `QPstartup.m` 预留位置加入：

```matlab
% add_mpaths_only([CPATH 'your_module/']);
```

## 在计算路径中挂接

物理主路径在 `packages/+qp/launcher.m`。若只需旁路模块，可在 `qp_driver.m` 调用前后，或 `service_driver.m` 中按 config 开关调用：

```matlab
if config.CONTROL.enable_your_module
  output.msg('v0s', '%s', 'Your module is enabled.');
  result = your_kernel(config);
end
```

新代码应通过 service manager 取数据（如 `system.get()`、`wave_functions.get()`）。

## 日志

使用 `service/+output`：

```matlab
cleanup = output.push('your_kernel');
output.msg('v1s', '%s', 'started');
```

建议层级：`v0` 主流程，`v1` 常规信息，`v2` 调试细节。
