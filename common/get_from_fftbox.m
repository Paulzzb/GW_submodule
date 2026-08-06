% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2025/06/26 ZZ

function data = get_from_fftbox(idxnz, fftbox, Nfft)
  % get_from_fftbox -> get data from fftbox

  ndata = length(idxnz);
  data = complex(zeros(ndata, 1), 0);
  if (any(size(fftbox) ~= Nfft))
    msg = ['Size of fftbox is not matched with Nfft!\nfftbox: ' ...
            num2str(size(fftbox)) ', Nfft: ' num2str(Nfft)];
    output.err(msg);
  end
  data(1:ndata) = fftbox(idxnz);
end

