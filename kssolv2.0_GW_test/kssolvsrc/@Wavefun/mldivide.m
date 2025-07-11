function Z = mldivide(X,R)
% WAVEFUN/MRDIVIDE Mrdivide function for wave function class
%    Z = MRDIVIDE(X,R) returns the wave function of X/R.
%
%    See also Wavefun.

Z = X.psi\R.psi;

end
