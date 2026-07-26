% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function coeff_coarse_fft_grid_build(id)
% Regular coarse FFT sublattice on a commensurate grid:
%   a = gcd(fftgrid_i), fftgrid_min = fftgrid_i/a, fftgrid_c = ratio*fftgrid_min
% with ratio minimal such that prod(fftgrid_c) >= nisdf target.

  fft_data = FFT.get();
  isdf_data = isdftest.get(id);

  fftgrid_i = int32(fft_data.fftgrid(:).');
  if numel(fftgrid_i) ~= 3
    error('isdf:coeff_coarse_fft_grid_build:BadFftgrid', 'fftgrid must have 3 components.');
  end

  Nmu_tmp = double(isdf_data.nisdf);

  nfft_d = double(fftgrid_i);
  total_fine = double(prod(nfft_d));
  if total_fine <= 0
    error('isdf:coeff_coarse_fft_grid_build:BadFine', 'invalid total fine-grid size.');
  end

  % Commensurate coarse lattice: fftgrid_c = ratio * fftgrid_min, with
  % fftgrid_min = fftgrid_i / a and a the common factor of fftgrid_i (gcd).
  % (Using lcm for a would give fftgrid_i/a not generally dividing the fine grid.)
  a = int32(local_gcd_int32(fftgrid_i));
  if a < 1
    error('isdf:coeff_coarse_fft_grid_build:BadGcd', 'invalid gcd of fftgrid_i.');
  end
  fftgrid_min = max(int32(1), fftgrid_i ./ a);

  ratio = int32(0);
  while prod(double(fftgrid_min) .* double(ratio)) <= Nmu_tmp
    ratio = ratio + 1;
  end
  fftgrid_c = fftgrid_min .* ratio;

  % for idim = 1:3
  %   if mod(fftgrid_i(idim), fftgrid_c(idim)) ~= 0
  %     error('isdf:coeff_coarse_fft_grid_build:NotCommensurate', ...
  %       'fftgrid_i(%d)=%d not divisible by fftgrid_c(%d)=%d (ratio=%d, a=%d).', ...
  %       idim, fftgrid_i(idim), idim, fftgrid_c(idim), ratio, a);
  %   end
  % end

  Nmu = int32(prod(double(fftgrid_c)));

  xq = (0:fftgrid_c(1) - 1);
  yq = (0:fftgrid_c(2) - 1);
  zq = (0:fftgrid_c(3) - 1);
  [Xq, Yq, Zq] = ndgrid(xq, yq, zq);
  R_coarse_RLU = double([Xq(:), Yq(:), Zq(:)]);
  R_rot_coarse = isdftest.coeff.isdf_coarse_R_rot(fftgrid_c);

  isdf_data.fftgrid_c = fftgrid_c;
  isdf_data.R_rot_coarse = R_rot_coarse;
  isdf_data.R_sampling_RLU = R_coarse_RLU;
  isdf_data.nisdf = Nmu;
  isdftest.save2mod(isdf_data, id);
end

function g = local_gcd_int32(v)
  g = int32(v(1));
  for k = 2:numel(v)
    g = int32(gcd(double(g), double(v(k))));
  end
end
