function Y = power(X,p)
% WAVEFUN/POWER Power function for wave function class
%    Y = POWER(X,p) returns the wave function of X.^p.
%
%    See also Wavefun.

Y = X;
Y.psi = X.psi.^p;

end
