function X = fft3(X)
% FFT3 Three-dimensional discrete Fourier Transform.
%     FFT3(X) returns the three-dimensional Fourier transform of tensor X.
%  
%     See also fft, fft2, fftn, ifft, ifft2, ifft3, ifftn.

X = fft(X,[],1);
X = fft(X,[],2);
X = fft(X,[],3);

end
