% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03 ZZ

function gen_coeff(cfg, id)
% ISDF index-controller function.
% This dispatcher selects an index-generation route and forwards control
% to the corresponding implementation function.

  if nargin < 1
    cfg = struct();
  end

  method = local_pick_method(cfg, id);

  switch lower(method)
    case 'default'
      isdf.coeff.gen_coeff_default(cfg, id);

    case 'qrcp'
      isdf.coeff.gen_coeff_qrcp(cfg, id);

    case 'kmeans'
      isdf.coeff.gen_coeff_kmeans(cfg, id);

    case 'coarse'
      isdf.coeff.gen_coeff_coarse(id);

    case 'pseudo'
      isdf.coeff.gen_coeff_pseudo(id);

    otherwise
      error(sprintf('ISDF.gen_coeff: unknown ISDF index method: %s', method));
  end
end

function method = local_pick_method(cfg, id)
  [field, default_method] = local_exxmethod_type_info(id);
  method = default_method;

  if ~isstruct(cfg)
    return;
  end

  if isfield(cfg, field) && ~isempty(cfg.(field))
    method = lower(strtrim(char(string(cfg.(field)))));
    local_warn_exxmethod_override(id, field, method, default_method);
    return;
  end

  if isfield(cfg, 'exxmethod') && ~isempty(cfg.exxmethod)
    method = lower(strtrim(char(string(cfg.exxmethod))));
    local_warn_exxmethod_override(id, 'exxmethod', method, default_method);
    return;
  end
end

function [field, default_method] = local_exxmethod_type_info(id)
  isdf_data = isdf.get(id);
  desc = lower(strtrim(char(string(isdf_data.desc))));

  switch desc
    case 'vc'
      field = 'exxmethod_type1';
      default_method = 'coarse';
    case 'vn'
      field = 'exxmethod_type2';
      default_method = 'pseudo';
    case 'nn'
      field = 'exxmethod_type3';
      default_method = 'coarse';
    otherwise
      field = 'exxmethod';
      default_method = 'coarse';
  end
end

function local_warn_exxmethod_override(id, field_name, method, default_method)
  if strcmp(method, default_method)
    return;
  end
  isdf_data = isdf.get(id);
  warning('isdf:gen_coeff:exxmethodOverride', ...
    ['ISDF id=%d (desc=%s): %s=''%s'' differs from default ''%s''.'], ...
    int32(id), char(string(isdf_data.desc)), field_name, method, default_method);
end
