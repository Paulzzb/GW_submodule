% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/24

function driver(data, config)
  % 
  [nk, nspin] = size(data.psig);
  nb = size(data.psig{1, 1}, 2);

  FFT_data = FFT.get();
  d_lat_data = lattice.manager('d_lat', 'get');
  DL_vol = d_lat_data.DL_vol;
  fftgrid = FFT_data.fftgrid;
  nc = prod(fftgrid);

  wf_data = wave_functions.base.wf_m(nc, nb, nk, nspin);

  for ispin = 1:nspin
    for ik = 1:nk
      data.psig{ik, ispin} = ...
      data.psig{ik, ispin} ./ sqrt(  sum( abs(data.psig{ik, ispin}).^2, 1 )  ) * sqrt(DL_vol);
    end
  end
  nknb = nk * nb;
  use_parfor = parallel.enabled();
  for ispin = 1:nspin
    c_spin_2d = complex(zeros(nc, nknb, 'double'));
    if use_parfor
      parfor ikib = 1:nknb
        ik = floor((ikib - 1) / nb) + 1;
        ib = mod(ikib - 1, nb) + 1;
        fftbox = put_into_fftbox(data.psig{ik, ispin}(:, ib), data.reciprocal_grid_info.idxnz{ik}, fftgrid);
        fftbox = nc ./ DL_vol * do_FFT(fftbox, fftgrid, 1);
        c_spin_2d(:, ikib) = double(fftbox(:));
      end
    else
      for ikib = 1:nknb
        ik = floor((ikib - 1) / nb) + 1;
        ib = mod(ikib - 1, nb) + 1;
        fftbox = put_into_fftbox(data.psig{ik, ispin}(:, ib), data.reciprocal_grid_info.idxnz{ik}, fftgrid);
        fftbox = nc ./ DL_vol * do_FFT(fftbox, fftgrid, 1);
        c_spin_2d(:, ikib) = double(fftbox(:));
      end
    end
    wf_data.c(:, :, :, ispin) = reshape(c_spin_2d, [nc, nb, nk]);
  end

  wave_functions.save2mod(wf_data);
end
