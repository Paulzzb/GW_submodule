% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function service_driver(data, config)
  % A temporary driver to call service modules for testing and demonstration.

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
    relay.save2db('test_relay_stage.mat', stage);
  catch ME
    warning(ME.identifier, 'Failed to save relay stage to DB: %s', ME.message);
  end
end
