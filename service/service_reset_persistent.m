% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19

function report = service_reset_persistent(opts)
  % Clear persistent caches across service modules.
  %
  %   service_reset_persistent(struct('keep_pool', true))
  %   skips parallel.free() so an existing parpool can be reused.

  keep_pool = false;
  if nargin >= 1 && ~isempty(opts) && isstruct(opts) ...
      && isfield(opts, 'keep_pool') && logical(opts.keep_pool)
    keep_pool = true;
  end

  cleanup_list = {
    % @() symmetry.free(), 'symmetry.free';
    % @() FFT.free(), 'FFT.free';
    % @() lattice.free(), 'lattice.free';
    % @() coulomb.free(), 'coulomb.free';
    % @() wave_functions.manager('free'), 'wave_functions.free';
    @() wave_functions.WF_apply_symm('reset'), 'wave_functions.WF_apply_symm(''reset'')';
    @() isdf.coeff.isdf_apply_symm_on_coarse('reset'), 'isdf.coeff.isdf_apply_symm_on_coarse(''reset'')';
    @() isdf.coeff.isdf_get_coeff('reset'), 'isdf.coeff.isdf_get_coeff(''reset'')';
    @() isdf.adaptive_double.isdf_schur_update('clear'), 'isdf.adaptive_double.isdf_schur_update(''clear'')';
    @() isdf.adaptive_double.adaptive_weight('clear'), 'isdf.adaptive_double.adaptive_weight(''clear'')';
    @() isdf.adaptive_single.isdf_schur_update('clear'), 'isdf.adaptive_single.isdf_schur_update(''clear'')';
    @() isdf.adaptive_single.adaptive_weight('clear'), 'isdf.adaptive_single.adaptive_weight(''clear'')';
    % get_u_xalpha('reset'): sampling2bundle/R_rot caches and u_xalpha_sampling_debug_on
    % (persistent cache of isdf.debug.on('u_xalpha_sampling'); see get_u_xalpha.m).
    @() isdf.get_u_xalpha('reset'), 'isdf.get_u_xalpha(''reset'')';
    @() isdf.get_rho_xalpha('reset'), 'isdf.get_rho_xalpha(''reset'')';
    @() isdf.debug.clear(), 'isdf.debug.clear';
    @() parallel.free(), 'parallel.free';
    @() timing.free(), 'timing.free';
    @() output.free(), 'output.free';
  };

  if keep_pool
    skip = {'parallel.free'};
    mask = true(size(cleanup_list, 1), 1);
    for k = 1:numel(skip)
      mask = mask & ~strcmp(cleanup_list(:, 2), skip{k});
    end
    cleanup_list = cleanup_list(mask, :);
  end

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
