function Y = ctranspose(X)
% WAVEFUN/CTRANSPOSE Ctranspose function for wave function class
%    Y = CTRANSPOSE(X) returns the conjugate transpose value of the wave
%    function.
%
%    See also Wavefun.

Y = X;
Y.trans = X.trans == 0;

end
