% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/01

function gen_coeff(cfg, id)
% ISDF index-controller function.
% This dispatcher selects an index-generation route and forwards control
% to the corresponding implementation function.
%
% Note:
% - This version only defines routing logic.
% - Concrete index-generation algorithms are intentionally left as stubs.

  if nargin < 1
    cfg = struct();
  end

  method = local_pick_method(cfg);

  switch lower(method)
    case 'default'
      isdftest.coeff.gen_coeff_default(cfg, id);

    case 'qrcp'
      isdftest.coeff.gen_coeff_qrcp(cfg, id);

    case 'kmeans'
      isdftest.coeff.gen_coeff_kmeans(cfg, id);

    case 'coarse'
      isdftest.coeff.gen_coeff_coarse(id);

    case 'pseudo'
      isdftest.coeff.gen_coeff_pseudo(id);

    otherwise
      error(sprintf('ISDF.gen_coeff: unknown ISDF index method: %s', method));
  end
end

function method = local_pick_method(cfg)
% Determine route method from config.

  method = 'coarse';

  if isstruct(cfg) && isfield(cfg, 'exxmethod') && ~isempty(cfg.exxmethod)
    method = char(string(cfg.exxmethod));
    return
  end

  if isstruct(cfg) && isfield(cfg, 'isdfoptions') && isstruct(cfg.isdfoptions)
    opts = cfg.isdfoptions;
    if isfield(opts, 'samp') && ~isempty(opts.samp)
      method = char(string(opts.samp));
      return
    end
  end
end
