% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/29 ZZ
function hVh = isdf_helper2hVh(helper, Dcoul, vol)

hVh = helper' * Dcoul * helper / vol;