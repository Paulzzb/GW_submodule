% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function rho_xalpha = get_rho_xalpha(id, param)
% GET_RHO_XALPHA  ISDF density on sampling points, then left-apply Lambda^{-ratio}V'.
%
% Computes c = phase * conj(u1) .* u2 as in +isdf, then returns
%   c_t = C_left(:,:,iqibz) * c, where C_left = Lambda^{-ratio}V' (padded).
% Downstream contraction is c_t' * tildeVq * c_t.

  persistent firsttime N_MAX inv_rot_index R_sampling fftgrid_c

  if nargin == 1 && (ischar(id) || (isstring(id) && isscalar(id)))
    cmd = lower(string(id));
    if cmd == "reset"
      firsttime = [];
      N_MAX = [];
      inv_rot_index = [];
      R_sampling = [];
      fftgrid_c = [];
      rho_xalpha = [];
      return;
    end
  end

  if isempty(N_MAX)
    N_MAX = isdf.isdf_nmax();
    firsttime = true(N_MAX, 1);
    symm_data = symmetry.get();
    inv_rot_index = symm_data.inv_rot_index;
    R_sampling = cell(N_MAX, 1);
    fftgrid_c = zeros(N_MAX, 3);
  end

  if firsttime(id)
    firsttime(id) = false;
    isdf_data = isdf.get(id);
    R_sampling{id} = isdf_data.bundle_struct.R_grid_bundle(isdf_data.bundle_struct.sampling2bundle, :);
    fft_data = FFT.manager('get');
    fftgrid_c(id, :) = fft_data.fftgrid;
  end

  if ~isfield(param, 'qs') || numel(param.qs) < 2
    error('isdf:get_rho_xalpha:qs', 'param.qs must be [iGo, iqibz, iqrot].');
  end

  r_lat_data = lattice.manager('r_lat', 'get');
  symm_data = symmetry.get();

  iGo = param.qs(1);
  iqibz = param.qs(2);
  iqrot = param.qs(3);

  u_xalpha1 = isdf.get_u_xalpha(id, param.is, iqrot);
  u_xalpha2 = isdf.get_u_xalpha(id, param.os, iqrot);
  c = conj(u_xalpha1) .* u_xalpha2;

  Go = double(r_lat_data.Ggrid_RLU(iGo, :));
  if norm(Go) >= 1e-6
    inviqrot = inv_rot_index(iqrot);
    invSqGo = double(Go * symm_data.rot_mtrx_RLU_G(:, :, inviqrot));
    phase_shift_coarse = exp(-2 * pi * 1i * (R_sampling{id} ./ double(fftgrid_c(id, :))) * invSqGo');
    c = c .* phase_shift_coarse;
  end

  isdf_data = isdf.get(id);
  has_cell_factors = isprop(isdf_data, 'CCHq_trunc_factors') && ...
                     numel(isdf_data.CCHq_trunc_factors) >= double(iqibz) && ...
                     ~isempty(isdf_data.CCHq_trunc_factors{double(iqibz)});
  if ~has_cell_factors
    error('isdf:get_rho_xalpha:MissingTruncFactors', ...
      ['CCHq_trunc_factors{%d} missing for id=%d; run isdf.gen_tildeVq first ', ...
       'and store V_trunc/Lambda_trunc per q.'], iqibz, int32(id));
  end
  fac = isdf_data.CCHq_trunc_factors{double(iqibz)};
  if ~isfield(fac, 'V_trunc') || ~isfield(fac, 'Lambda_trunc')
    error('isdf:get_rho_xalpha:BadFactors', ...
      'CCHq_trunc_factors{%d} missing V_trunc/Lambda_trunc.', iqibz);
  end
  isdf.prod_C_inv_t('set_factors', fac.V_trunc, fac.Lambda_trunc, isdf_data.svd_ratio);
  rho_xalpha = isdf.prod_C_inv_t('prod', 'l', c(:));
end
