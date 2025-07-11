function X = compress(X,idxnz)
% WAVEFUN/COMPRESS Compress function for wave function class
%    X = COMPRESS(X,idxnz) returns the compact wave function with the
%    non-zero index idxnz.
%
%    See also Wavefun.

X.idxnz = idxnz;
X.psi = X.psi(idxnz,:);
X.iscompact = 1;

end
