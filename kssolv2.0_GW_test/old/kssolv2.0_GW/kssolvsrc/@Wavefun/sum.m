function beta = sum(X,varargin)
% WAVEFUN/SUM Sum function for wave function class
%    beta = SUM(X) returns the sum of the wave function.
%
%    beta = SUM(X,dim) returns the sum of the wave function along the dim
%    dimension.
%
%    See also Wavefun.

beta = sum(X.psi,varargin{:});

end
