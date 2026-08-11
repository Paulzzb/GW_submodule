% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function showtag(flag)
%SHOWTAG  Enable/disable [module] prefixes on messages.
%
%   output.showtag(true)
%   output.showtag(false)

  if nargin < 1 || ~(islogical(flag) || isnumeric(flag))
    error('output:showtag:BadFlag', 'showtag expects logical true/false.');
  end
  s = state_('get');
  s.showtag = logical(flag);
  state_('set', s);
end
