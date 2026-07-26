% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06

function tf = enabled()
  try
    obj = parallel.get();
    tf = obj.enabled;
  catch
    tf = license('test', 'Distrib_Computing_Toolbox');
  end
end
