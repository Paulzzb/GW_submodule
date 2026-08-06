% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ

function summary(stage, config, varargin)
%SUMMARY  Staged QP summary (pre- and post-calculation).
%
%   qp.summary(0, config)              % before compute
%   qp.summary(1, config, E, elapsed)  % after compute

  cleanup = output.push('+qp/summary.m'); %#ok<NASGU>

  if nargin < 2
    output.err('qp.summary(stage, config, ...) requires stage and config.');
  end

  stage = double(stage);
  switch stage
    case 0
      stage0(config);
    case 1
      if numel(varargin) < 2
        output.err('qp.summary(1, config, E, elapsed) requires E and elapsed.');
      end
      stage1(config, varargin{1}, varargin{2});
    otherwise
      output.err('qp.summary stage must be 0 (pre) or 1 (post), got %g.', stage);
  end
end

function stage0(config)
  output.msg('v0s', 'Eqp = Ex + Esx_x + Ecoh');

  freq_dep = config.FREQUENCY.frequency_dependence;
  switch freq_dep
    case -2
      route = 'gw.x_Gamma + gw.cohsex_Gamma';
    case -1
      route = 'gw.x + gw.cohsex_multi_k';
    case 2
      route = 'gw.fullfreq_cd_res_Gamma + gw.fullfreq_cd_int_Gamma';
    otherwise
      route = sprintf('unknown (frequency_dependence=%g)', freq_dep);
  end
  output.msg('v0s', 'Route: frequency_dependence=%g → %s', freq_dep, route);

  output.warn('Multi-spin is not supported yet.');
  if freq_dep == -2 || freq_dep == 2
    output.warn('Only Gamma point is supported.');
  elseif freq_dep == -1
    output.warn('Multi-k path is a prototype (ikibz/ispin fixed in places).');
  end

  if ~isfield(config, 'ISDF') || ~logical(config.ISDF.isisdf)
    output.msg('v0s', 'ISDF: off (dense G-space)');
  else
    ex_label = "vn";
    if isfield(config, 'COHSEX') && isfield(config.COHSEX, 'ex_use_which_isdf') ...
        && ~isempty(config.COHSEX.ex_use_which_isdf)
      ex_label = lower(string(config.COHSEX.ex_use_which_isdf));
    end

    try
      [id_vc, id_vn, id_nn] = isdf.resolve_ids(config);
      output.msg('v0s', 'ISDF: %s; %s; %s (ex_use_which_isdf=%s)', ...
        local_id_str('vc', id_vc), local_id_str('vn', id_vn), ...
        local_id_str('nn', id_nn), char(ex_label));
    catch ME
      output.warn('ISDF slot resolve failed: %s', ME.message);
    end
  end

  if isfield(config, 'freqinfo') && ~isempty(config.freqinfo) ...
      && isfield(config.freqinfo, 'grid_real')
    nreal = numel(config.freqinfo.grid_real);
    nimag = numel(config.freqinfo.grid_imag);
    res_m = config.FREQUENCY.cd_residual_method;
    int_m = config.FREQUENCY.cd_integration_method;
    output.msg('v0s', ...
      'Frequency: %d real / %d imag (methods %d / %d)', ...
      nreal, nimag, res_m, int_m);
  end
end

function s = local_id_str(label, id)
  if isempty(id)
    s = sprintf('%s=<missing>', label);
  else
    s = sprintf('%s=%d', label, double(id));
  end
end

function stage1(~, E, elapsed)
  nb = 0;
  nk = 0;
  if isfield(E, 'Eqp') && ~isempty(E.Eqp)
    nb = size(E.Eqp, 1);
    nk = size(E.Eqp, 2);
  end
  output.msg('v0s', 'Finished in %.2f s (nb=%d, nk=%d).', ...
    double(elapsed), nb, nk);

  if isfield(E, 'fout') && ~isempty(E.fout)
    output.msg('v0s', 'Saved table to: %s', char(string(E.fout)));
  end
end
