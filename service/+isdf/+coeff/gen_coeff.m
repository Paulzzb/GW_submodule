% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function idx_mu = gen_coeff(cfg, id)
% ISDF index-controller function.
% This dispatcher selects an index-generation route and returns fine-grid
% linear indices. Slot fill / bundle init is done by the caller via
% isdf.rsymm.init_from_indices(id, idx_mu).
%
%   idx_mu = isdf.coeff.gen_coeff(cfg, id)
%
% Also writes interp_scheme on the ISDF slot (method metadata).
% For method 'coarse', the legacy gen_coeff_coarse path still fills the slot
% itself and this function returns empty indices.
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
      output.err('ISDF.gen_coeff: unknown ISDF index type: %s', desc);
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
      idx_mu = isdf.coeff.gen_coeff_default(cfg, id);
      scheme = "default";

    case 'qrcp'
      idx_mu = isdf.coeff.gen_coeff_qrcp(cfg, id);
      scheme = "qrcp";

    case 'kmeans'
      idx_mu = isdf.coeff.gen_coeff_kmeans(cfg, id);
      scheme = "kmeans";

    case 'coarse'
      % Legacy path: fills the slot itself; no indices returned for
      % init_from_indices. Not part of the indices-only +coeff redesign yet.
      isdf.coeff.gen_coeff_coarse(id);
      idx_mu = int32([]);
      return;

    case 'pseudo'
      idx_mu = isdf.coeff.gen_coeff_pseudo(id);
      scheme = "pseudo";

    otherwise
      output.err('ISDF.gen_coeff: unknown ISDF index method: %s', method);
  end

  idx_mu = int32(idx_mu(:));
  isdf_data = isdf.get(id);
  isdf_data.interp_scheme = scheme;
  isdf.save2mod(isdf_data, id);
end
