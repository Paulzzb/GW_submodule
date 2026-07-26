% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function id = isdftest_add(desc)
% Allocate first free cell in the ISDF pool (N_MAX = 10), set desc, select its id.
% Returns int32 id (analogous spirit to Yambo FFT_add).

  if nargin < 1 || isempty(desc)
    desc = "";
  end
  id = isdftest.manager('add', desc);
end
