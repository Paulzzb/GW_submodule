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
    % get_u_xalpha caches R_rot / sampling2bundle per id in persistent state; after pool
    % teardown the same numeric id may refer to a new ISDF slot with a fresh bundle.
    isdf.get_u_xalpha('reset');
  else
    info = isdf.manager('free', id);
  end
end
