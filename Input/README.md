### what is doing
构建如下的工作流
[inputfile] → read_input_param() → struct(config) → convert_to_GWOptions(config) → GWOptions 实例（支持 display & 快捷访问）
仍在构建中...

## Current framework
GWOptionsFramework/
Input/                           # Input parsing and validation layer
├── 
│   ├── read_input_param.m           # Read and parse the input file into a config struct
│   ├── allowed_param_list.m         # Define valid blocks and parameter names
│   ├── default_param_values.m       # Provide default values for each parameter
│   ├── set_default_values.m         # Fill missing fields in config using defaults
│   └── ...
│
├── @GWOptions/                      # Core class representing parameter configuration
│   ├── GWOptions.m                  # Class definition with constructor
│   ├── fromConfig.m                 # Convert struct (config) to GWOptions object
│   ├── display.m                    # Customized display of essential parameters
│   └── ...                          # Additional methods like validation or shortcuts
│
├── examples/                        # Example input files (optional)
│   └── input_example.in
│
└── README.md                        # Project documentation (this file)
## Purpose
现在的GW module的输入都是直接提供一个struct options_GW，然后修改。
但是实际上这样即使我一个源代码开发者都搞的乱七八糟
现在我期待，重新弄一套规范化输入。
- [ ] 将GW module里面的关键参数分成几大类，
每一大类对应着某一个关键功能。
- [ ] 提供一个稳定的、可扩展的参数读取器。
- [ ] 用markdown书写一个面相用户的参数手册。
- [ ] 测试

Transform the Input of other codes,
such as QE for the ground-state calculation, or BGW for GW calculation,
to the inputs form the code can accept.
## Input description
