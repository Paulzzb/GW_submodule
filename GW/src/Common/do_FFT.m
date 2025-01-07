function fftbox = do_FFT(fftbox, Nfft, sign)
    % The 'fftbox' parameter in MATLAB is already complex, so no separate real and imaginary parts are needed.
    % MATLAB's FFT functions directly work on complex data.

    if sign == 1
        fftbox = ifftn(fftbox, Nfft);  % Inverse FFT (backward transform)
    elseif sign == -1
        fftbox = fftn(fftbox, Nfft);   % Forward FFT
    else
        error('sign is not 1 or -1 in do_FFT');
    end
end

