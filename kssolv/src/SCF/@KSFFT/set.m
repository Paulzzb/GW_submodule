function F = set(F,varargin)
% KSFFT/SET Set function for KSFFT class
%    F = SET(F,str1,field1,str2,field2,...) returns a KSFFT class of
%    the given fields with respect to the name strings.
%
%    See also KSFFT.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Hargin = varargin;
while (length(Hargin)>=2),
   attr_name = Hargin{1};
   value     = Hargin{2};
   Hargin    = Hargin(3:end);
   switch attr_name
     case 'n1'
        F.n1 = value;
     case 'n2'
        F.n2 = value; 
     case 'n3'           
        F.n3 = value;    
     case 'vol'
        F.vol = vol;
     case 'idxnz'
        F.idxnz = value;
     case 'forward'
        F.forward= value;
     case 'inverse'
        F.inverse = value;
     otherwise
        fprintf('Attribute name: %s\n', attr_name);
        error('cannot be found');
   end;
end;
