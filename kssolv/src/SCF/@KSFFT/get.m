function val = get(F,attr_name)
% KSFFT/GET Obtain attribute of KSFFT class
%    val = GET(F,attr_name) returns attribute in F
%
%    See also KSFFT.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch attr_name
  case 'n1'
    val = F.n1;
  case 'n2'
    val = F.n2;
  case 'n3'
    val = F.n3;
  case 'vol'
    val = F.vol;
  case 'idxnz'
    val = F.idxnz;
  case 'forward'
    val = F.forward;
  case 'inverse'
    val = F.inverse;
  otherwise
    error('invalid attribute requested for KSFFT');
end;
