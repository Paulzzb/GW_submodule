function Xabs = abs(X)
% WAVEFUN/ABS Abs function for wave function class
%    Xabs = ABS(X) returns the absolute value of the wave function.
%
%    See also Wavefun.

Xabs = X;
Xabs.psi = abs(X.psi);

end
