% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ

function opts = cauchy_opts(config)
%CAUCHY_OPTS  Thin Cauchy options from config.ISDF.iscauchy.
%
%   opts = isdf.cauchy_opts(config)
%   fields: isCauchy (logical), froErr, MaxIter (fixed defaults)

  opts = struct('isCauchy', false, 'froErr', 1e-6, 'MaxIter', 10);
  if nargin < 1 || isempty(config) || ~isfield(config, 'ISDF') ...
      || ~isstruct(config.ISDF)
    return
  end

  if isfield(config.ISDF, 'iscauchy') && ~isempty(config.ISDF.iscauchy)
    opts.isCauchy = logical(config.ISDF.iscauchy);
  elseif isfield(config.ISDF, 'isCauchy') && ~isempty(config.ISDF.isCauchy)
    opts.isCauchy = logical(config.ISDF.isCauchy);
  end
end
