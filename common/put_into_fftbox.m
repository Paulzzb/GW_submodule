% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2025/06/13 ZZ

function fftbox = put_into_fftbox(data, idxnz, Nfft)
% put_into_fftbox -> put data into fftbox

ndata = length(data);
fftbox = complex(zeros(Nfft(1), Nfft(2), Nfft(3)), 0);  
fftbox(idxnz) = data(1:ndata);

end % EOF

