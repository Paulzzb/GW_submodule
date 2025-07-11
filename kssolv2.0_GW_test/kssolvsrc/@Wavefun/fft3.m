function X = fft3(X)
% WAVEFUN/FFT3 FFT3 function for wave function class
%    Xfft = FFT3(X) returns the fast Fourier transform of the wave
%    function.
%
%    See also Wavefun.

if X.iscompact
    xpsi = zeros(X.n1*X.n2*X.n3,ncols(X));
    xpsi(X.idxnz,:) = X.psi;
    xpsi = reshape(xpsi,X.n1,X.n2,X.n3,ncols(X));
    X.idxnz = [];
else
    xpsi = reshape(X.psi,X.n1,X.n2,X.n3,ncols(X));
end

xpsi = fft3(xpsi);
X.psi = reshape(xpsi,X.n1*X.n2*X.n3,ncols(X));
X.iscompact = 0;

end
