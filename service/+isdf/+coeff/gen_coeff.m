% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/05 ZZ

function gen_coeff(cfg, id)
% ISDF index-controller function.
% This dispatcher selects an index-generation route and forwards control
% to the corresponding implementation function.
%
% Default method for all types is 'pseudo'. Enum strings on cfg are assumed
% normalized by set_default_param_value (trim + lowercase).

  if nargin < 1 || isempty(cfg)
    cfg = struct();
  end

  isdf_data = isdf.get(id);
  desc = lower(strtrim(char(string(isdf_data.desc))));

  switch desc
    case 'vc'
      field = 'exxmethod_type1';
    case 'vn'
      field = 'exxmethod_type2';
    case 'nn'
      field = 'exxmethod_type3';
    otherwise
      output.error('ISDF.gen_coeff: unknown ISDF index type: %s', desc);
  end
  default_method = 'pseudo';

  if isfield(cfg, field) && ~isempty(cfg.(field))
    method = char(string(cfg.(field)));
  else
    method = default_method;
  end

  if ~strcmp(method, default_method)
    output.warn('ISDF id=%d (desc=%s): %s=''%s'' differs from default ''%s''.', ...
      int32(id), char(string(isdf_data.desc)), field, method, default_method);
  end

  switch method
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
      error('ISDF.gen_coeff: unknown ISDF index method: %s', method);
  end
end
