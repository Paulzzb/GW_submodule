% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function cfg = get_relay_config()
% Central relay configuration: module name => list of persistent variables to track.
% Modify this function alone to add/remove tracked modules and variables.
%
% Intentionally NOT staged (rebuilt every service_driver call from current config):
%   lattice, coulomb, ISDF

  cfg = struct();
  cfg.version = 2;

  % Define modules and their persistent variables
  % Format: cfg.modules.<module_name>.<var_name> = struct('type', '<manager_call_signature>')

  % ========== symmetry ==========
  cfg.modules.symmetry = struct();
  cfg.modules.symmetry.main = struct();
  cfg.modules.symmetry.main.type = 'double';
  cfg.modules.symmetry.main.call_get = @() symmetry.manager('get');
  cfg.modules.symmetry.main.call_set = @(val) symmetry.manager('save2mod', val);

  % ========== system ==========
  cfg.modules.system = struct();
  cfg.modules.system.main = struct();
  cfg.modules.system.main.type = 'double';
  cfg.modules.system.main.call_get = @() system.manager('get');
  cfg.modules.system.main.call_set = @(val) system.manager('save2mod', val);

  % ========== FFT ==========
  cfg.modules.FFT = struct();
  cfg.modules.FFT.main = struct();
  cfg.modules.FFT.main.type = 'double';
  cfg.modules.FFT.main.call_get = @() FFT.manager('get');
  cfg.modules.FFT.main.call_set = @(val) FFT.manager('save2mod', val);

  % ========== pair_symmetry ==========
  cfg.modules.pair_symmetry = struct();
  cfg.modules.pair_symmetry.main = struct();
  cfg.modules.pair_symmetry.main.type = 'double';
  cfg.modules.pair_symmetry.main.call_get = @() pair_symmetry.manager('get');
  cfg.modules.pair_symmetry.main.call_set = @(val) pair_symmetry.manager('save2mod', val);

  % ========== wave_functions ==========
  cfg.modules.wave_functions = struct();
  cfg.modules.wave_functions.main = struct();
  cfg.modules.wave_functions.main.type = 'double';
  cfg.modules.wave_functions.main.call_get = @() wave_functions.manager('get');
  cfg.modules.wave_functions.main.call_set = @(val) wave_functions.manager('save2mod', val);

  % ========== timing ==========
  cfg.modules.timing = struct();
  cfg.modules.timing.main = struct();
  cfg.modules.timing.main.type = 'double';
  cfg.modules.timing.main.call_get = @() timing.manager('get');
  cfg.modules.timing.main.call_set = @(val) timing.manager('save2mod', val);

end
