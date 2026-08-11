% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03

function params = adaptive_param(desc, cfg_isdf)
%ADAPTIVE_PARAM  Map config.ISDF (+ desc) into adaptive run params.
%
%   params = isdf.adaptive.adaptive_param(desc, cfg_isdf)
%
% Defaults are owned by input/default_param_values.m (+ set_default_param_value).
% This function does not invent fallback numbers; cfg_isdf must already be filled.

  if nargin < 2 || isempty(cfg_isdf) || ~isstruct(cfg_isdf)
    error('adaptive_param:cfg', ...
      ['cfg_isdf (config.ISDF) is required. ', ...
       'Initialize once via default_param_values / set_default_param_value, then pass config.ISDF.']);
  end

  desc = lower(strtrim(char(string(desc))));
  if strcmp(desc, 'vc')
    suffix = 'type1';
  elseif strcmp(desc, 'vn')
    suffix = 'type2';
  elseif strcmp(desc, 'nn')
    suffix = 'type3';
  else
    error('adaptive_param:desc', ...
      'Unsupported ISDF desc ''%s'' (expect vc/vn/nn).', desc);
  end

  params = struct();
  params.threshold = local_require_positive(cfg_isdf, ['adaptive_threshold_' suffix]);
  params.num_add = int32(local_require_integer(cfg_isdf, ['adaptive_num_add_' suffix]));
  params.candidate_ratio = local_require_positive(cfg_isdf, ['adaptive_candidate_ratio_' suffix]);
  params.isdf_ratio = local_require_positive(cfg_isdf, ['isdf_ratio_' suffix]);
  params.max_cond = local_require_positive(cfg_isdf, ['adaptive_max_cond_' suffix]);
  params.use_cond_guard = local_require_logical(cfg_isdf, 'adaptive_use_cond_guard');
  params.weight_batch_size = local_resolve_batch_size(cfg_isdf);
  params.source = "config";
end

function val = local_require_positive(cfg, field_name)
  if ~isfield(cfg, field_name)
    error('adaptive_param:field', 'Missing config.ISDF.%s', field_name);
  end
  v = double(cfg.(field_name));
  if ~(isfinite(v) && v > 0)
    error('adaptive_param:field', ...
      'config.ISDF.%s must be a positive finite value (got %g).', field_name, v);
  end
  val = v;
end

function val = local_require_integer(cfg, field_name)
  if ~isfield(cfg, field_name)
    error('adaptive_param:field', 'Missing config.ISDF.%s', field_name);
  end
  v = double(cfg.(field_name));
  if ~(isfinite(v) && v >= 1)
    error('adaptive_param:field', ...
      'config.ISDF.%s must be an integer >= 1 (got %g).', field_name, v);
  end
  val = max(1, round(v));
end

function val = local_require_logical(cfg, field_name)
  if ~isfield(cfg, field_name)
    error('adaptive_param:field', 'Missing config.ISDF.%s', field_name);
  end
  v = cfg.(field_name);
  if islogical(v)
    val = logical(v);
    return;
  end
  if isnumeric(v) && isfinite(v) && isscalar(v)
    val = logical(v ~= 0);
    return;
  end
  if ischar(v) || isstring(v)
    s = lower(strtrim(char(string(v))));
    if any(strcmp(s, {'true', '.true.', '1', 'yes', 'on'}))
      val = true;
      return;
    elseif any(strcmp(s, {'false', '.false.', '0', 'no', 'off'}))
      val = false;
      return;
    end
  end
  error('adaptive_param:field', ...
    'config.ISDF.%s must be a logical (or 0/1).', field_name);
end

function val = local_resolve_batch_size(cfg)
  if ~isfield(cfg, 'adaptive_batch_size')
    error('adaptive_param:field', 'Missing config.ISDF.adaptive_batch_size');
  end
  v = double(cfg.adaptive_batch_size);
  if ~(isfinite(v) && v >= 1)
    error('adaptive_param:field', ...
      'config.ISDF.adaptive_batch_size must be an integer >= 1 (got %g).', v);
  end
  val = max(1, round(v));
end
