% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11 ZZ

function attach_gamma_trunc_factors(id, inv_ratio)
%ATTACH_GAMMA_TRUNC_FACTORS  Identity SVD factors after gen_tildeVq_Gamma.
%
% Gamma builders skip the SVD path in gen_tildeVq; gw.x_Gamma / gen_Kq_Gamma
% still expect CCHq_trunc_factors. This sets full-rank identity factors.

  isdf_data = isdf.get(id);
  if isempty(isdf_data.tildeVq)
    error('isdf:attach_gamma_trunc_factors:NoTildeVq', ...
      'ISDF id=%d has empty tildeVq; run gen_tildeVq_Gamma first.', int32(id));
  end
  Nmu = size(isdf_data.tildeVq, 1);
  nibz = size(isdf_data.tildeVq, 3);
  if nargin < 2 || isempty(inv_ratio)
    inv_ratio = 0.75;
  end
  if ~isfinite(inv_ratio) || inv_ratio <= 0
    inv_ratio = 0.75;
  end
  isdf_data.svd_ratio = inv_ratio;
  isdf_data.CCHq_trunc_factors = cell(1, double(nibz));
  for iq = 1:double(nibz)
    fac = struct();
    fac.V_trunc = eye(Nmu);
    fac.Lambda_trunc = ones(Nmu, 1);
    fac.N_keep = int32(Nmu);
    isdf_data.CCHq_trunc_factors{iq} = fac;
  end
  isdf.save2mod(isdf_data, id);
end
