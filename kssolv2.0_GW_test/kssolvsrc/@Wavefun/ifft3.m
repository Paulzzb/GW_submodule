function Xfft = ifft3(Xfft)
% WAVEFUN/IFFT3 IFFT3 function for wave function class
%    X = IFFT3(Xfft) returns the inverse fast Fourier transform of the wave
%    function.
%
%    See also Wavefun.

if Xfft.iscompact
    fpsi = zeros(Xfft.n1*Xfft.n2*Xfft.n3,ncols(Xfft));
    fpsi(Xfft.idxnz,:) = Xfft.psi;
    fpsi = reshape(fpsi,Xfft.n1,Xfft.n2,Xfft.n3,ncols(Xfft));
    Xfft.idxnz = [];
else
    fpsi = reshape(Xfft.psi,Xfft.n1,Xfft.n2,Xfft.n3,ncols(Xfft));
end

fpsi = ifft3(fpsi);
Xfft.psi = reshape(fpsi,Xfft.n1*Xfft.n2*Xfft.n3,ncols(Xfft));
Xfft.iscompact = 0;

end
