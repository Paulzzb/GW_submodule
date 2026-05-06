% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function coeff_coarse_fft_grid_build(id)
% Regular coarse FFT sublattice: prod(fftgrid_c) > isdf_ratio * nb (wave_functions.nb).

  fft_data = FFT.get();
  isdf_data = isdf.get(id);

  fftgrid_i = int32(fft_data.fftgrid(:).');
  if numel(fftgrid_i) ~= 3
    error('isdf:coeff_coarse_fft_grid_build:BadFftgrid', 'fftgrid must have 3 components.');
  end

  Nmu_tmp = single(isdf_data.nisdf);

  nfft_d = single(fftgrid_i);
  total_fine = single(prod(nfft_d));
  if total_fine <= 0
    error('isdf:coeff_coarse_fft_grid_build:BadFine', 'invalid total fine-grid size.');
  end

  scale = (max(Nmu_tmp, 1) / total_fine)^(1 / 3);
  fftgrid_c = int32(max(1, ceil(scale .* nfft_d)));
  while prod(single(fftgrid_c)) <= Nmu_tmp
    cand = zeros(1, 3);
    for idim = 1:3
      d_try = single(fftgrid_c);
      d_try(idim) = d_try(idim) + 1;
      cand(idim) = prod(d_try);
    end
    [~, pick] = min(cand);
    fftgrid_c(pick) = fftgrid_c(pick) + int32(1);
  end

  Nmu = int32(prod(single(fftgrid_c)));

  xq = (0:fftgrid_c(1) - 1);
  yq = (0:fftgrid_c(2) - 1);
  zq = (0:fftgrid_c(3) - 1);
  [Xq, Yq, Zq] = ndgrid(xq, yq, zq);
  R_coarse_RLU = single([Xq(:), Yq(:), Zq(:)]);
  R_rot_coarse = isdf.coeff.isdf_coarse_R_rot(fftgrid_c);

  isdf_data.fftgrid_c = fftgrid_c;
  isdf_data.R_rot_coarse = R_rot_coarse;
  isdf_data.R_sampling_RLU = R_coarse_RLU;
  isdf_data.nisdf = Nmu;
  isdf.save2mod(isdf_data, id);
end
