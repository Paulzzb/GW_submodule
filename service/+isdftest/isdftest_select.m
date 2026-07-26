% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function isdftest_select(id)
% Set current ISDF pool id for isdftest.get / isdftest.save2mod (no isolation; global current pointer).

  isdftest.manager('select', id);
end
