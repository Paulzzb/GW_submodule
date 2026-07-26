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

  def = filename_map();
  dir = config.CONTROL.storage_dir;
  fNameGWinput = fullfile(dir, def.GWinput);
  use_cache = isfile(fNameGWinput);

  if use_cache
    fprintf('input_driver: loading cached SAVE data from %s\n', dir);
    GWgroundstate = load(fNameGWinput, 'GWgroundstate').GWgroundstate;
    cached = load(fullfile(dir, def.config), 'GWoptions', 'config');
    GWoptions = cached.GWoptions;
    config = cached.config;
    data = load(fullfile(dir, def.data), 'data').data;
    % Fill defaults for fields added after config.mat was saved.
    config = set_default_param_value(config, data);
  else
    % Step 4: Load groundstate info, transform them into a uniform format
    %         struct 'data'.
    dirin = config.CONTROL.groundstate_dir;
    typein = config.CONTROL.groundstate_type;
    data = load_groundstate_info(dirin, typein, config);

    % Step 3: Set default values in 'config', error if there exists invalid values.
    config = set_default_param_value(config, data);

    % Step 4: Construct GWinfo (Basically, the groundstate data) and GWOptions seperately
    GWgroundstate = construct_GWinfo(data, config);
    GWoptions = construct_GWOptions(data, config);

    % Full-frequency (contour deformation): frequency grids for gw_fullfreq_cd_* / qp_cohsex.
    if config.FREQUENCY.frequency_dependence == 2
      config = generate_frequency(GWgroundstate, config);
    end

    % Step 5: save data to files
    fNamedata = fullfile(dir, def.data);
    fName2 = fullfile(dir, def.config);
    if ~exist(dir, 'dir')
      mkdir(dir);
    end

    config.ISDFCauchy = GWoptions.ISDFCauchy;
    save(fNamedata, 'data', '-v7.3', '-nocompression');
    save(fNameGWinput, 'GWgroundstate', '-v7.3', '-nocompression');
    save(fName2, 'GWoptions', 'config');
  end
  %
  % Step 8: ( testing )
  % use structure in service/ to construct 
  service_driver(data, config);
  % Step 7: display input and groundstate information
  display_input_summary(GWgroundstate, GWoptions, config)

end % function