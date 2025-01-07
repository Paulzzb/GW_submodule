% function fftbox = put_into_fftbox(ndata, data, gvec)
function fftbox = put_into_fftbox(data, idxnz, Nfft)
% Inputs
% ndata: number of data
% data
% glist: size at least (ndata, 3);
% gindex: to put it inside.
% Nfft: size of fftbox
% Output:
% fftbox: with proper number inside.

% if nargin < 6
% 	error('method is not given in put_into_fftbox.')
% end

im                 = sqrt(-1);
% CZERO              = 0.0 + im * 0.0;

ndata = length(data);

fftbox = complex(zeros(Nfft(1), Nfft(2), Nfft(3)), 0);  

% for j = 1:ndata
%     bidx = gvec_to_fft_index(glist(gindex(j), :), Nfft);
%     fftbox(bidx(1), bidx(2), bidx(3)) = data(j);
% end
fftbox(idxnz) = data(1:ndata);


return
end
