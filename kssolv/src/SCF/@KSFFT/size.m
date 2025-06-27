function [nrow, ncol] = size(F)
% KSFFT/SIZE The information of mesh size for KSFFT
%    [nrow, ncol] = SIZE(F) returns the mesh size for 
%    FFT and the number of G grids
%
%    See also KSFFT.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nrow = 0;
ncol = 0;
if (F.inverse)
   nrow = F.n1*F.n2*F.n3;
   ncol = length(F.idxnz);
elseif (F.forward)
   nrow = length(F.idxnz); 
   ncol = F.n1*F.n2*F.n3;
end
