% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [Nmu, ind_mu] = gen_coarse_Rgrid()
% Coarse subset of fine FFT R-grid sites (gcd / step rule).
% Nmu_tmp = isdf_ratio * nrep; pick smallest regular sublattice with Nmu > Nmu_tmp.

  fft_data = FFT.get();
  isdf_data = isdf.get();
  pair_data = pair_symmetry.get();

  fftgrid = int32(fft_data.fftgrid(:).');
  if numel(fftgrid) ~= 3
    error('isdf:coeff_coarse_rgrid_indices:BadFftgrid', 'fftgrid must have 3 components.');
  end

  nrep = double(pair_data.nrep);
  Nmu_tmp = double(isdf_data.isdf_ratio) * nrep;

  g = int32(gcd(gcd(fftgrid(1), fftgrid(2)), fftgrid(3)));
  if g <= 0
    error('isdf:coeff_coarse_rgrid_indices:BadGcd', 'invalid gcd(fftgrid) = %d.', g);
  end

  base = int32(fftgrid ./ g);

  m_list = isdf.coeff.divisors_int32(g);
  m_list = sort(m_list, 'ascend');

  chosen_m = int32(-1);
  chosen_Nmu = -1;
  for i = 1:numel(m_list)
    m = m_list(i);
    coarse_dims = int32(base .* m);
    nmu_candidate = double(prod(double(coarse_dims)));
    if nmu_candidate > Nmu_tmp
      chosen_m = m;
      chosen_Nmu = nmu_candidate;
      break;
    end
  end

  if chosen_m < 1
    chosen_m = m_list(end);
    chosen_Nmu = double(prod(double(base .* chosen_m)));
  end

  step = int32(g / chosen_m);
  if step < 1
    error('isdf:coeff_coarse_rgrid_indices:BadStep', 'invalid step size computed.');
  end

  rgrid = fft_data.Rgrid_RLU;
  if size(rgrid, 2) ~= 3
    error('isdf:coeff_coarse_rgrid_indices:BadRgrid', 'Rgrid_RLU must be nr x 3.');
  end

  rgrid_i = int32(round(double(rgrid)));
  mask = mod(rgrid_i(:, 1), step) == 0 & ...
         mod(rgrid_i(:, 2), step) == 0 & ...
         mod(rgrid_i(:, 3), step) == 0;

  ind_mu = int32(find(mask));
  Nmu = int32(numel(ind_mu));

  if double(Nmu) ~= chosen_Nmu
    error('isdf:coeff_coarse_rgrid_indices:NmuMismatch', ...
      'selected Nmu mismatch (%d vs %d).', Nmu, int32(chosen_Nmu));
  end
end
