function w = ncols(BX)
% BLOCHWAVEFUN/NCOLS Numbers of columns in the Bloch wave function class
%    w = NCOLS(BX) returns the numbers of columns in the Bloch wave
%    function.
%
%    See also BlochWavefun.

w = cellfun(@(X)ncols(X),BX.wavefuncell);

end
