% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function input_driver(inputfile)
%INPUT_DRIVER  Prepare SAVE data and service-layer objects for a GW run.
%
%   Cache gate is SAVE/relay_stage.mat (expensive). Config is always rebuilt
%   from the namelist + defaults (cheap) so CUTOFFS/ISDF edits take effect
%   without deleting the stage; lattice/coulomb/ISDF are rebuilt in
%   service_driver from that fresh config.

  % Step 1: Read and parse namelist -> fresh config
  config = read_input_param(inputfile);

  % Step 2: Validate required fields
  validate_required_params(config);

  def = filename_map();
  dir = config.CONTROL.storage_dir;
  fNameStage = fullfile(dir, def.stage);
  fNamedata = fullfile(dir, def.data);
  fNameconfig = fullfile(dir, def.config);
  use_stage = isfile(fNameStage) && isfile(fNamedata);

  if use_stage
    fprintf('input_driver: stage cache hit (%s); reloading data, rebuilding config\n', ...
      fNameStage);
    data = load(fNamedata, 'data').data;
  else
    if isfile(fNameStage) && ~isfile(fNamedata)
      warning('input_driver:StageWithoutData', ...
        ['Found %s but missing %s; ignoring stage and rebuilding from groundstate.'], ...
        fNameStage, fNamedata);
    end
    dirin = config.CONTROL.groundstate_dir;
    typein = config.CONTROL.groundstate_type;
    data = load_groundstate_info(dirin, typein, config);
  end

  % Always (re)derive config from namelist + data defaults — never reuse config.mat.
  config = set_default_param_value(config, data);

  % Open r-<prefix>.log before anything that may output.msg('r', ...), including
  % generate_frequency. Otherwise a leftover FID from a previous run appends
  % those lines to the end of the old report, which is then rotated to r-*_NN.
  output.free();
  display_input_summary(config);

  if config.FREQUENCY.frequency_dependence == 2
    config = generate_frequency(config);
  end

  if ~exist(dir, 'dir')
    mkdir(dir);
  end
  config.ISDFCauchy = setISDFCauchy(data, config);
  save(fNamedata, 'data', '-v7.3', '-nocompression');
  save(fNameconfig, 'config', '-v7.3', '-nocompression');

  service_driver(data, config);

end % function
