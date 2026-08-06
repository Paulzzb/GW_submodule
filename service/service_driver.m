% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06

function service_driver(data, config)
  % A temporary driver to call service modules for testing and demonstration.

  def = filename_map();
  storage_dir = config.CONTROL.storage_dir;
  stage_path = fullfile(storage_dir, def.stage);

  if isfile(stage_path)
    fprintf('service_driver: loading cached relay stage from %s\n', stage_path);
    relay.stage_from_db(stage_path);
    relay.restore();
    if config.ISDF.isisdf
      % Relay only snapshots isdf.manager('get') (one current slot), not the full
      % vc/vn/nn pool built by isdf.driver. Without rebuilding, qp.launcher cannot
      % resolve ISDF ids after a cache hit.
      isdf.driver(data, config);
    end
    return
  end

  isdf.debug.init_from_config(config);
  parallel.driver(data, config);

  system.driver(data, config);
  symmetry.driver(data, config);
  FFT.driver(data, config);
  lattice.driver(data, config);
  coulomb.driver(data, config);
  % pair_symmetry.driver(data, config);
  wave_functions.driver(data, config);
  if config.ISDF.isisdf
    isdf.driver(data, config);
  end

  stage = relay.collect();
  try
    if ~exist(storage_dir, 'dir')
      mkdir(storage_dir);
    end
    relay.save2db(stage_path, stage);
  catch ME
    warning(ME.identifier, 'Failed to save relay stage to DB: %s', ME.message);
  end
end
