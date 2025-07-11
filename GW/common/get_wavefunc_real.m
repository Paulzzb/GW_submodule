function psir = get_wavefunc_real(psig, Ggrid_info)

[nkpts, nspin] = size(psig);


fftgrid = Ggrid_info.fftgrid;
nfftgrid = prod(fftgrid);
vol = Ggrid_info.vol;
nb = size(psig{1,1}, 2);

for ispin = 1:nspin
  for ikpt = 1:nkpts
    for ib = 1:nb
      psig{ikpt, ispin}(:, ib) = psig{ikpt, ispin}(:, ib) / norm(psig{ikpt, ispin}(:, ib));
    end
    psig{ikpt, ispin} = psig{ikpt, ispin} * sqrt(vol);
  end
end


psir = zeros(nfftgrid, nb, nkpts, nspin);
for ispin = 1:nspin
  for ikpt = 1:nkpts
    for iband = 1:nb
      fftbox = put_into_fftbox(psig{ikpt, ispin}(:, iband), Ggrid_info.idxnz{ikpt}, fftgrid);
      fftbox = nfftgrid ./ vol * do_FFT(fftbox, fftgrid, 1);
      psir(:, iband, ikpt, ispin) = fftbox(:);
    end
  end
end


end % EOF