% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

function obj = get(id)
% isdf.get()  current id;  isdf.get(id)  read that id without changing current.

  if nargin < 1
    obj = isdf.manager('get');
  else
    obj = isdf.manager('get', id);
  end
end
