% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20 ZZ

function input_driver(inputfile)
% This is a driver function, given an inputfile, it will prepare and save all
% required data for GW module calculation in storage_dir.
% The code does the following step:
% 1. Read the input file, generate a struct 'config'.
% 2. Validate required parameters.
% 3. If cached SAVE data exists in storage_dir, load it and fill new defaults;
%    otherwise:
%    a. Load groundstate data from groundstate_dir into a uniform struct 'data'.
%    b. Set default values in 'config' from defaults and groundstate data.
%    c. Construct GWgroundstate and GWOptions from 'data' and 'config'.
%    d. For full-frequency (contour deformation), generate frequency grids.
%    e. Save data, GWgroundstate, GWOptions and config to storage_dir.
% 4. Build service-layer objects via service_driver.
% 5. Display input and groundstate summary.
  

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
    % GWgroundstate = load(fNameGWinput, 'GWgroundstate').GWgroundstate;
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
    % GWgroundstate = construct_GWinfo(data, config);
    GWoptions = construct_GWOptions(data, config);

    % Full-frequency (contour deformation): frequency grids for gw_fullfreq_cd_* / qp_cohsex.
    if config.FREQUENCY.frequency_dependence == 2
      config = generate_frequency(config);
    end

    % Step 5: save data to files
    fNamedata = fullfile(dir, def.data);
    if ~exist(dir, 'dir')
      mkdir(dir);
    end

    config.ISDFCauchy = setISDFCauchy(data, config);
    save(fNamedata, 'data', '-v7.3', '-nocompression');
  end
  %
  % Step 8: ( testing )
  % use structure in service/ to construct 
  service_driver(data, config);
  % Step 7: display input and groundstate information
  display_input_summary(config)

end % function