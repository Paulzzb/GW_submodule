function Xconj = conj(X)
% WAVEFUN/CONJ Conj function for wave function class
%    Xconj = CONJ(X) returns the conjugate value of the wave function.
%
%    See also Wavefun.

Xconj = X;
Xconj.psi = conj(X.psi);

end
