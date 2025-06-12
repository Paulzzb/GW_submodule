function aqstemp = mtxel_sigma(nn, GWinfo, sum_range)
% mtxel_sigma -> Calculate <nn | exp(iGr) | mm>, where mm in sum_range
%     use GWinfo.psir, GWinfo.gvec, and GWinfo.vol
% if sum_range is not specified, sum_range = size(GWinfo.psir, 2)
% 

if nargin < 3
  sum_range = size(GWinfo.psir, 2);
end

vol = GWinfo.vol;
gvec = GWinfo.gvec;

aqstemp = complex(0.0, 0.0) * zeros(ng, length(sum_range));


for ind = 1:length(sum_range)
  fftbox = conj(GWinfo.psir(:, nn)) .* GWinfo.psir(:, sum_range(ind));
  fftbox = reshape(fftbox, GWinfo.gvec.fftgrid);
  fftbox = vol * do_FFT(fftbox, gvec.fftgrid, 1);
  aqstemp(:, ind) = get_from_fftbox(gvec.idxnz, fftbox, gvec.fftgrid);
end

end % EOF 
