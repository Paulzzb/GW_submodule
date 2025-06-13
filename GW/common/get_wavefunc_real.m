function psir = get_wavefunc_real(psig, Ggrid_info)

fftgrid = Ggrid_info.fftgrid;
nfftgrid = prod(fftgrid);
vol = Ggrid_info.vol;
nb = size(psig, 2);


for ib = 1:nb
  psig(:, ib) = psig(:, ib) / norm(psig(:, ib));
end
psig = psig * sqrt(vol);

psir = zeros(nfftgrid, nb);
for iband = 1:nb
  fftbox = put_into_fftbox(psig(:, iband), Ggrid_info.idxnz, fftgrid);
  fftbox = nfftgrid ./ vol * do_FFT(fftbox, fftgrid, 1);
  psir(:, iband) = fftbox(:);
end

end % EOF