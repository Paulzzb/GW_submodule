% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function cfg = get_relay_config()
% Central relay configuration: module name => list of persistent variables to track.
% Modify this function alone to add/remove tracked modules and variables.

  cfg = struct();
  cfg.version = 1;

  % Define modules and their persistent variables
  % Format: cfg.modules.<module_name>.<var_name> = struct('type', '<manager_call_signature>')

  % ========== symmetry ==========
  cfg.modules.symmetry = struct();
  cfg.modules.symmetry.main = struct();
  cfg.modules.symmetry.main.type = 'single';
  cfg.modules.symmetry.main.call_get = @() symmetry.manager('get');
  cfg.modules.symmetry.main.call_set = @(val) symmetry.manager('save2mod', val);

  % ========== system ==========
  cfg.modules.system = struct();
  cfg.modules.system.main = struct();
  cfg.modules.system.main.type = 'single';
  cfg.modules.system.main.call_get = @() system.manager('get');
  cfg.modules.system.main.call_set = @(val) system.manager('save2mod', val);

  % ========== FFT ==========
  cfg.modules.FFT = struct();
  cfg.modules.FFT.main = struct();
  cfg.modules.FFT.main.type = 'single';
  cfg.modules.FFT.main.call_get = @() FFT.manager('get');
  cfg.modules.FFT.main.call_set = @(val) FFT.manager('save2mod', val);

  % ========== lattice ==========
  cfg.modules.lattice = struct();
  cfg.modules.lattice.r_lat = struct();
  cfg.modules.lattice.r_lat.type = 'keyed';
  cfg.modules.lattice.r_lat.key = 'r_lat';
  cfg.modules.lattice.r_lat.call_get = @() lattice.manager('r_lat', 'get');
  cfg.modules.lattice.r_lat.call_set = @(val) lattice.manager('r_lat', 'save2mod', val);

  cfg.modules.lattice.d_lat = struct();
  cfg.modules.lattice.d_lat.type = 'keyed';
  cfg.modules.lattice.d_lat.key = 'd_lat';
  cfg.modules.lattice.d_lat.call_get = @() lattice.manager('d_lat', 'get');
  cfg.modules.lattice.d_lat.call_set = @(val) lattice.manager('d_lat', 'save2mod', val);

  cfg.modules.lattice.k = struct();
  cfg.modules.lattice.k.type = 'keyed';
  cfg.modules.lattice.k.key = 'k';
  cfg.modules.lattice.k.call_get = @() lattice.manager('k', 'get');
  cfg.modules.lattice.k.call_set = @(val) lattice.manager('k', 'save2mod', val);

  cfg.modules.lattice.q = struct();
  cfg.modules.lattice.q.type = 'keyed';
  cfg.modules.lattice.q.key = 'q';
  cfg.modules.lattice.q.call_get = @() lattice.manager('q', 'get');
  cfg.modules.lattice.q.call_set = @(val) lattice.manager('q', 'save2mod', val);

  % ========== coulomb ==========
  cfg.modules.coulomb = struct();
  cfg.modules.coulomb.main = struct();
  cfg.modules.coulomb.main.type = 'single';
  cfg.modules.coulomb.main.call_get = @() coulomb.manager('get');
  cfg.modules.coulomb.main.call_set = @(val) coulomb.manager('save2mod', val);

  % ========== pair_symmetry ==========
  cfg.modules.pair_symmetry = struct();
  cfg.modules.pair_symmetry.main = struct();
  cfg.modules.pair_symmetry.main.type = 'single';
  cfg.modules.pair_symmetry.main.call_get = @() pair_symmetry.manager('get');
  cfg.modules.pair_symmetry.main.call_set = @(val) pair_symmetry.manager('save2mod', val);

  % ========== wave_functions ==========
  cfg.modules.wave_functions = struct();
  cfg.modules.wave_functions.main = struct();
  cfg.modules.wave_functions.main.type = 'single';
  cfg.modules.wave_functions.main.call_get = @() wave_functions.manager('get');
  cfg.modules.wave_functions.main.call_set = @(val) wave_functions.manager('save2mod', val);

  % ========== ISDF ==========
  cfg.modules.ISDF = struct();
  cfg.modules.ISDF.main = struct();
  cfg.modules.ISDF.main.type = 'single';
  cfg.modules.ISDF.main.call_get = @() isdf.manager('get');
  cfg.modules.ISDF.main.call_set = @(val) isdf.manager('save2mod', val);

  % ========== timing ==========
  cfg.modules.timing = struct();
  cfg.modules.timing.main = struct();
  cfg.modules.timing.main.type = 'single';
  cfg.modules.timing.main.call_get = @() timing.manager('get');
  cfg.modules.timing.main.call_set = @(val) timing.manager('save2mod', val);


end
