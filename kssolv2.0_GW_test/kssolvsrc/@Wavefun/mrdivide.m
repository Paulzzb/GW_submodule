function Z = mrdivide(X,R)
% WAVEFUN/MRDIVIDE Mrdivide function for wave function class
%    Z = MRDIVIDE(X,R) returns the wave function of X/R.
%
%    See also Wavefun.

Z = X;
Z.psi = X.psi/R;

end
