function Y = uminus(X)
% WAVEFUN/UMINUS Uminus function for wave function class
%    Y = UMINUS(X) returns a wave function as -X.
%
%    See also Wavefun.

Y = X;
Y.psi = -X.psi;

end
