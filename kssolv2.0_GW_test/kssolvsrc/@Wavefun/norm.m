function beta = norm(X,varargin)
% WAVEFUN/NORM Norm function for wave function class
%    beta = NORM(X,opt) returns the norm of the wave functions.
%
%    See also Wavefun.

beta = norm(X.psi,varargin{:});

end
