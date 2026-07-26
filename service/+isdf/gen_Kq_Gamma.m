% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/08 ZZ
function Kq_ISDF = gen_Kq_Gamma(id_vc, iqibz, omega)
% This function only consider single-k point sampling.
% With the stable adaptiveISDF

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

if nargin < 3
  omega = complex(0);
end

eta = 0.0;
warning('Currently, only consider omega = 0 case.');
nspin = 1; ispin = 1;
warning('Multi-spin is not supported yet.');

system_data = system.get();
k_data = lattice.manager('k', 'get');
r_lat_data = lattice.manager('r_lat', 'get');
% 
nkbz = k_data.nbz;
ev = system_data.Eo;
% ev = system_data.Eo * ry2ev;

vc_data = isdf.get(id_vc);
if isempty(vc_data.nrange1) || isempty(vc_data.nrange2)
  error('gen_Kq:nrange', ...
    'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', ...
    int32(id_vc));
end
nrangev = double(vc_data.nrange1);
nrangec = double(vc_data.nrange2);
Nisdf_vc = vc_data.nisdf;
s2b = vc_data.bundle_struct.sampling2bundle;

Nisdf_o = Nisdf_vc;
R_sampling = vc_data.bundle_struct.R_grid_bundle(vc_data.bundle_struct.sampling2bundle, :);
fft_data = FFT.manager('get');
fftgrid = fft_data.fftgrid;
symm_data = symmetry.get();
inv_rot_index = symm_data.inv_rot_index;


iqrot = 1;
ikibz = 1;
chiq_ISDF = zeros(Nisdf_o, Nisdf_o);
nrangec_row = reshape(nrangec, 1, []);
f_c = system_data.f(nrangec_row, 1, ispin);
e_c = ev(nrangec_row, 1, ispin);
% Now we always assume that occupation numbers are either 0 or 1.
for iv = nrangev
  f_ik = system_data.f(iv, ikibz, ispin);
  e_ik = ev(iv, ikibz, ispin);
  occ = -f_c + f_ik;
  if eta == 0.0
    den = omega - e_ik + e_c;
  else
    den = omega - e_ik + e_c + 1i * eta * sign(e_ik - e_c);
  end
  valid = (abs(occ) >= 1e-5) & (abs(den) >= 1e-12);
  if ~any(valid)
    continue;
  end

  rho_blk = zeros(Nisdf_o, nnz(valid));
  coeff = occ(valid) ./ den(valid);
  %
  is = [iv, 1, 1, ispin];
  u_xalpha_is = isdf.get_u_xalpha(id_vc, is, iqrot);
  %
  u_xalpha_os = vc_data.bundle_struct.WF_bundle(s2b, nrangec_row, 1, ispin);
  rho_blk = conj(u_xalpha_is) .* u_xalpha_os;
  %
  rho_weighted = rho_blk .* reshape(coeff, 1, []);
  chiq_ISDF = chiq_ISDF + (rho_weighted * rho_blk');
end
% when occupation number are either 0 or 1, the extra factor 2.0 is needed (ij-pair and ji-pair).
chiq_ISDF = 2.0 * chiq_ISDF;
% when spin == 1, an extra factor 2.0 is needed (up-down and down-up pairs).
if nspin == 1
  chiq_ISDF = 2 * chiq_ISDF;
end
% symmetrize the matrix if real(omega) = 0.
if abs(real(omega)) < 1e-12
  chiq_ISDF = (chiq_ISDF + chiq_ISDF') / 2;
end

% Beware that we need to multiply \Lambda^{inv_ratio-1}_t * V_t
inv_ratio = vc_data.svd_ratio;
fac = vc_data.CCHq_trunc_factors{double(iqibz)};
V_trunc = fac.V_trunc;
Lambda_trunc = fac.Lambda_trunc;
N_keep = fac.N_keep;

chiq_ISDF = chiq_ISDF * V_trunc;
chiq_ISDF = V_trunc' * chiq_ISDF;
chiq_ISDF = Lambda_trunc.^(inv_ratio-1) .* (chiq_ISDF) .* Lambda_trunc.^(inv_ratio-1).';
chiq_ISDF = 0.5 * chiq_ISDF' + 0.5 * chiq_ISDF;

% Then K
% K_q(\mu, \nu) = \Lambda_q^{-1}(\mu, \nu) - \tildeV_q(\mu, \nu) * inv_ratio
Kq_ISDF = inv(chiq_ISDF) + vc_data.tildeVq(1:N_keep, 1:N_keep, iqibz);

end