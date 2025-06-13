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

  
  % Step 4: Load groundstate info, transform them into a uniform format
  %         struct 'data'.
  dirin = config.CONTROL.groundstate_dir;
  typein = config.CONTROL.groundstate_type;
  data = load_groundstate_info(dirin, typein);

  % Step 3: Set default values in 'config', error if there exists invalid values.
  config = set_default_param_value(config, data);


  % Step 4: Construct GWinfo (Basically, the groundstate data) and GWOptions seperately
  GWgroundstate = construct_GWinfo(data, config);
  GWoptions = construct_GWOptions(data, config);

  % Step 6: save data to files 
  %% [optionsGW, GWinfor] = construct_GW_objects(config, GWinfor);% 
  dir = config.CONTROL.storage_dir;
  def = filename_map();
  fName = fullfile(dir, def.GWinput);
  if ~exist(dir, 'dir')
    mkdir(dir);
  end
  
  config.ISDFCauchy = GWoptions.ISDFCauchy;
  save(fName, 'GWgroundstate', 'GWoptions', 'config');

  % Step 7: display input and groundstate information
  display_input_summary(GWgroundstate, GWoptions, config)

end % function