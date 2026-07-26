% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

function info = save2mod(data, id)
% isdf.save2mod(data)  current id;  isdf.save2mod(data, id)  explicit id.

  if nargin < 2
    info = isdf.manager('save2mod', data);
  else
    info = isdf.manager('save2mod', data, id);
  end
end
