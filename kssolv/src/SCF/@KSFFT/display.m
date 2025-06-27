function display(F)
% KSFFT/DISPLAY Prints the information of the KSFFT class
%    DISPLAY(F) shows the mesh size in real space, the
%    grid number in reciprocal space and whether the FFT
%    is forward or backward.
%
%    See also KSFFT.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

fprintf('KSFFT configuration: \n');
fprintf('n1 = %d, n2 = %d, n3 = %d\n', F.n1, F.n2, F.n3);
fprintf('ng = %d\n', length(F.idxnz));
if (F.forward) 
   fprintf('forward\n');
else
   fprintf('inverse\n');
end;
