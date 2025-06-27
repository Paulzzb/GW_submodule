function val = get(a,attr_name)
% ATOM/GET Obtain attribute of Atom class
% 	val = GET(a,attr_name) returns attribute in a
%
%   See also Atom.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch attr_name
  case 'symbol'
    val = mol.symbol;
  case 'anum'
    val = mol.anum;
  case 'amass'
    val = mol.amass;
  case 'venum'
    val = mol.venum;
  otherwise 
    error('invalid attribute requested for Atom');
end
