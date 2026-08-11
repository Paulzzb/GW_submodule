% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function idnew = launcher(id, cfg_isdf)
%LAUNCHER  Route adaptive ISDF to single or double backend by slot desc / threshold.
%
%   idnew = isdf.adaptive.launcher(id, cfg_isdf)
%
% Routing (cutoff = 1e-6):
%   vc           -> adaptive_single
%   vn / nn      -> adaptive_single if adaptive_threshold > cutoff, else adaptive_double
%
% Threshold fields: adaptive_threshold_type1/2/3 for vc/vn/nn.

  if nargin < 2 || isempty(cfg_isdf) || ~isstruct(cfg_isdf)
    error('adaptive:launcher:cfg', ...
      ['cfg_isdf (config.ISDF) is required. ', ...
       'Initialize once via default_param_values / set_default_param_value.']);
  end

  cutoff = constant_map().ADAPTIVE_THRESHOLD_CUTOFF;
  isdf_data = isdf.get(id);
  desc = lower(strtrim(char(string(isdf_data.desc))));

  switch desc
    case 'vc'
      idnew = isdf.adaptive_single.adaptiveisdf(id, cfg_isdf);

    case 'vn'
      thr = cfg_isdf.adaptive_threshold_type2;
      if thr <= cutoff
        idnew = isdf.adaptive_double.adaptiveisdf(id, cfg_isdf);
      else
        idnew = isdf.adaptive_single.adaptiveisdf(id, cfg_isdf);
      end

    case 'nn'
      thr = cfg_isdf.adaptive_threshold_type3;
      if thr <= cutoff
        idnew = isdf.adaptive_double.adaptiveisdf(id, cfg_isdf);
      else
        idnew = isdf.adaptive_single.adaptiveisdf(id, cfg_isdf);
      end

    otherwise
      error('adaptive:launcher:desc', ...
        'Unsupported ISDF desc ''%s'' (expect vc/vn/nn).', desc);
  end
end

