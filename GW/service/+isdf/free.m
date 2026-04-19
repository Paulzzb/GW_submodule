% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

function info = free(id)
% isdf.free()  clear entire pool;  isdf.free(id)  clear one id.

  if nargin < 1
    info = isdf.manager('free');
  else
    info = isdf.manager('free', id);
  end
end
