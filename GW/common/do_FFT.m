function fftbox = do_FFT(fftbox, Nfft, sign)
% do_FFT - Perform FFT on a box

  if sign == 1
    fftbox = ifftn(fftbox, Nfft);  % Inverse FFT (backward transform)
  elseif sign == -1
    fftbox = fftn(fftbox, Nfft);   % Forward FFT
  else
    msg = 'sign is not 1 or -1 in do_FFT';
    GWerror(msg);
  end
end

