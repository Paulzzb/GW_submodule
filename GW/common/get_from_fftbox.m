function data = get_from_fftbox(idxnz, fftbox, Nfft)
  % get_from_fftbox -> get data from fftbox

  ndata = length(idxnz);
  data = complex(zeros(ndata, 1), 0);
  if (any(size(fftbox) ~= Nfft))
    msg = ['Size of fftbox is not matched with Nfft!\nfftbox: ' ...
            num2str(size(fftbox)) ', Nfft: ' num2str(Nfft)];
    GWerror(msg);
  end
  data(1:ndata) = fftbox(idxnz);
end

