function input_driver(inputfile)
% This is a driver function, given an inputfile, it will save all required data
% for GW module calculation in storage_dir.
% The code does the following step:
% 1. Read the input file, generate a struct 'config'.
% 2. Check if necessary but missed parameters.
% 3. Load groundstate data from groundstate_dir.
% 4. Set default values with default_param_map, then with groundstate data.
% 5. Construct @GWinfo and @GWOptions based on groundstate data and 'config'.
  
  % Step 1: Read and parse
  config = read_input_param(inputfile);
  
  % Step 2: Validate required fields
  validate_required_params(config);
  
  % Step 3: Load groundstate info
  dirin = config.CONTROL.groundstate_dir;
  typein = config.CONTROL.groundstate_type;
  GWinfor = load_groundstate_info(dirin, typein);
  
  % Step 4: Fill in defaults (standard + based on GWinfor)
  config = complete_config_with_defaults(config, GWinfor);
  
  % Step 5: Build GWOptions and GWinfo
  [optionsGW, GWinfor] = construct_GW_objects(config, GWinfor);% 

end % function