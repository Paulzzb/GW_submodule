% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/24

function report = service_reset_persistent()
  % Clear persistent caches across service modules.

  cleanup_list = {
    % @() symmetry.free(), 'symmetry.free';
    % @() FFT.free(), 'FFT.free';
    % @() lattice.free(), 'lattice.free';
    % @() coulomb.free(), 'coulomb.free';
    % @() wave_functions.manager('free'), 'wave_functions.free';
    @() wave_functions.WF_apply_symm('reset'), 'wave_functions.WF_apply_symm(''reset'')';
    @() isdf.isdf_apply_symm_on_coarse('reset'), 'isdf.isdf_apply_symm_on_coarse(''reset'')';
    @() isdf_get_coeff('reset'), 'isdf_get_coeff(''reset'')';
    @() isdf_schur_update('clear'), 'isdf_schur_update(''clear'')';
    @() adaptive_weight('clear'), 'adaptive_weight(''clear'')';
    @() timing.free(), 'timing.free';
  };



  report.ok = true;
  report.cleared = {};
  report.errors = {};

  for i = 1:size(cleanup_list, 1)
    fn = cleanup_list{i, 1};
    name = cleanup_list{i, 2};

    try
      fn();
      report.cleared{end + 1} = name; %#ok<AGROW>
    catch ME
      report.ok = false;
      report.errors{end + 1} = sprintf('%s: %s', name, ME.message); %#ok<AGROW>
    end
  end
end
