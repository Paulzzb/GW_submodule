function data = get_from_fftbox(idxnz, fftbox, Nfft)
% Inputs
% ndata: number of data
% glist: size at least (ndata, 3);
% fftbox: Where to get data.
% gindex: to put it inside.
% Nfft: size of fftbox
%
% Output:
% data: size (ndata, 1)
  ndata = length(idxnz);
  data = complex(zeros(ndata, 1), 0);
  if (any(size(fftbox) ~= Nfft))
    error('Size of fftbox is not matched with Nfft!') 
  end
  % for j = 1:ndata
  %     bidx = gvec_to_fft_index(glist(gindex(j), :), Nfft);
  %     data(j) = fftbox(bidx(1), bidx(2), bidx(3));
  % end
  data(1:ndata) = fftbox(idxnz);
end

