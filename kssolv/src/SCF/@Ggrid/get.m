function val = get(grid,attr_name)
% GGRID/GET Obtain attribute of Ggrid class
%    val = GET(grid,attr_name) returns attribute in grid
%
%    See also Ggrid.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch attr_name
  case 'ecut'
    val = grid.ecut;
  case 'ng'
    val = grid.ng;
  case 'idxnz'
    val = grid.idxnz;
  case 'gkk'
    val = grid.gkk;
  case 'gkx'
    val = grid.gkx;
  case 'gky'
    val = grid.gky;
  case 'gkz'
    val = grid.gkz;
  case 'kkxyz'
    val = grid.kkxyz;
  otherwise
    error('Invalid attribute requested for Ggrid');
end;
