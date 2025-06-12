function fftbox = put_into_fftbox(data, idxnz, Nfft)
% put_into_fftbox -> put data into fftbox

ndata = length(data);
fftbox = complex(zeros(Nfft(1), Nfft(2), Nfft(3)), 0);  
fftbox(idxnz) = data(1:ndata);

end % EOF

