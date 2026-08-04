% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06

function info = free()
  if license('test', 'Distrib_Computing_Toolbox')
    try
      pool = gcp('nocreate');
      if ~isempty(pool)
        delete(pool);
      end
    catch ME
      warning('parallel:free:pool', ...
        'Failed to close parallel pool cleanly: %s', ME.message);
    end
  end
  info = parallel.manager('free');
end
