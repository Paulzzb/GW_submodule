% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19

function report = service_reset_persistent()
  % Clear persistent caches across service modules.

  cleanup_list = {
    % @() symmetry.free(), 'symmetry.free';
    % @() FFT.free(), 'FFT.free';
    % @() lattice.free(), 'lattice.free';
    % @() coulomb.free(), 'coulomb.free';
    % @() wave_functions.manager('free'), 'wave_functions.free';
    @() wave_functions.WF_apply_symm('reset'), 'wave_functions.WF_apply_symm(''reset'')';
    @() isdf.coeff.isdf_apply_symm_on_coarse('reset'), 'isdf.coeff.isdf_apply_symm_on_coarse(''reset'')';
    @() isdf.coeff.isdf_get_coeff('reset'), 'isdf.coeff.isdf_get_coeff(''reset'')';
    @() isdf.adaptive.isdf_schur_update('clear'), 'isdf.adaptive.isdf_schur_update(''clear'')';
    @() isdf.adaptive.adaptive_weight('clear'), 'isdf.adaptive.adaptive_weight(''clear'')';
    @() isdf.get_u_xalpha('reset'), 'isdf.get_u_xalpha(''reset'')';
    @() isdf.get_rho_xalpha('reset'), 'isdf.get_rho_xalpha(''reset'')';
    @() isdf.debug.clear(), 'isdf.debug.clear';
    @() parallel.free(), 'parallel.free';
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
