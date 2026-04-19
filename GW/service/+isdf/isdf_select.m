% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function isdf_select(id)
% Set current ISDF pool id for isdf.get / isdf.save2mod (no isolation; global current pointer).

  isdf.manager('select', id);
end
