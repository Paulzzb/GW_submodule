function Y = ctranspose(X)
% KSFFT/CTRANSPOSE Ctranspose function for KSFFT class
%    Y = CTRANSPOSE(X) returns a forward Fourier transform if X
%    is a backward Fourier transform, returns a backward Fourier 
%    transform if X is a forward Fourier transform
%
%    See also KSFFT.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin == 1)
   if (X.forward) 
      X.forward = 0;
      X.inverse = 1;
   elseif (X.inverse)
      X.inverse = 0;
      X.forward = 1; 
   else
      error('not sure which direction the original transform goes');
   end;
   Y = X;
else
   error('KSFFT: invalid syntax')
end;
