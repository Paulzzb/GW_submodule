% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function wf_on_coarse = coeff_coarse_wf_extract(fftgrid_c, R_coarse_RLU)
% Wavefunctions on coarse RLU sites: fine real-space WF -> FFT -> G_table -> dense map to coarse RLU.

  fft_data = FFT.get();
  fftgrid_i = int32(fft_data.fftgrid(:).');
  fft_sz = double(fftgrid_i);
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');

  DL_vol = d_lat_data.DL_vol;
  Nmu = int32(prod(double(fftgrid_c(:)')));

  scale_vec = double(fftgrid_i) ./ double(fftgrid_c(:)');
  R_scaled = double(R_coarse_RLU) .* scale_vec ./ double(fftgrid_i);
  phase_arg = 2 * pi * (R_scaled * double(r_lat_data.Ggrid_RLU)');
  scal = 1 ./ (prod(fft_sz) * double(DL_vol));
  iFFTmat = exp(1i * phase_arg) * scal;

  nibz = int32(k_data.nibz);
  nb = int32(wf_data.nb);
  nspin = int32(wf_data.nspin);

  wf_on_coarse = zeros(Nmu, nb, nibz, nspin);

  for ik_ibz = 1:nibz
    ik_rot = 1;
    for ispin = 1:nspin
      for ib = 1:nb
        isc = int32([ib, ik_ibz, ik_rot, ispin]);
        wf_b = wave_functions.WF_apply_symm(isc);
        fftbox = double(reshape(wf_b, fft_sz));
        fftbox = do_FFT(fftbox, fft_sz, -1) * double(DL_vol);
        psig = fftbox(fft_data.G_table(:, 1));
        wf_b_coarse = (iFFTmat * psig).';
        wf_on_coarse(:, ib, ik_ibz, ispin) = wf_b_coarse(:);
      end
    end
  end

end
