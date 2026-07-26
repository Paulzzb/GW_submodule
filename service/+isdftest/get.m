% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

function obj = get(id)
% isdftest.get()  current id;  isdftest.get(id)  read that id without changing current.

  if nargin < 1
    obj = isdftest.manager('get');
  else
    obj = isdftest.manager('get', id);
  end
end
