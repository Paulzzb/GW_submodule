function display(grid)
% GGRID/DISPLAY Print the number of truncated G grids 
%
%    DISPLAY(grid) shows the ecut and the grid number in
%    reciprocal space
%
%    See also Ggrid.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

fprintf('Energy cutoff used: %11.3e\n', grid.ecut);
ng = grid.ng;
if ( iscell(ng) )
   numkpts = length(ng);
   for kpt = 1:numkpts
      fprintf('kpoint %d:\n', kpt);
      fprintf('   number of nonzeros within ecut: %d\n', grid.ng{kpt});
   end
else
   % no a solid, no k point
   fprintf('number of nonzeros within ecut: %d\n', grid.ng);
   fprintf('          gx             gy              gz             |g|^2\n');
   for i = 1:grid.ng
      fprintf('%15.6e %15.6e %15.6e %15.6e\n', ...
              grid.gkx(i), grid.gky(i),grid.gkz(i), grid.gkk(i));
   end
end;
