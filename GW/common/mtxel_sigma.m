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
ng = gvec.ng;
idxnz = gvec.idxnz;
fftgrid = gvec.fftgrid;
nfft = prod(fftgrid);

aqstemp = complex(0.0, 0.0) * zeros(ng, length(sum_range));


for ind = 1:length(sum_range)
  if ~isempty(GWinfo.psir)
    fftbox1 = conj(GWinfo.psir(:, nn)) .* GWinfo.psir(:, sum_range(ind));
    fftbox1 = reshape(fftbox1, fftgrid);
  else
    fftbox1 = put_into_fftbox(GWinfo.psig(:, nn), idxnz, fftgrid) * sqrt(vol);
    fftbox1 = nfft / vol * do_FFT(fftbox1, fftgrid, 1);
    fftbox1 = conj(fftbox1);
    fftbox2 = put_into_fftbox(GWinfo.psig(:, sum_range(ind)), idxnz, fftgrid) * sqrt(vol);
    fftbox2 = nfft / vol * do_FFT(fftbox2, fftgrid, 1);
    fftbox1 = fftbox1.*fftbox2;
  end
  fftbox1 = vol * do_FFT(fftbox1, fftgrid, 1);
  aqstemp(:, ind) = get_from_fftbox(idxnz, fftbox1, fftgrid);
end

end % EOF 
