% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function level = verbose(new_level)
%VERBOSE  Get or set the global verbosity level (0/1/2).
%
%   n = output.verbose()
%   output.verbose(1)

  s = state_('get');
  if nargin < 1
    level = s.verbose;
    return
  end
  if ~isnumeric(new_level) || ~isscalar(new_level)
    error('output:verbose:BadLevel', 'verbose expects a scalar integer.');
  end
  s.verbose = max(0, min(2, round(double(new_level))));
  state_('set', s);
  level = s.verbose;
end
