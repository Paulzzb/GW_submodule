function Kq_ISDF = gen_Kq_Gamma(id_vc, iqibz, omega)

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

if nargin < 3
  omega = complex(0);
end

%   \chi_q(\mu, \nu)
% = 2 * \sum_{k\in\bz} \sum_{ij} * (f_{j\k-\q}-f_{i\k})
%     * \frac{1}{omega-Eo_{ik}+Eo_{j\k-\q} + i*eta*sign(Eo_{ik}-Eo_{j\k-\q})}
%     * \rho_{ij}(k, q, r_mu) \conj{\rho_{ij}(k, q, r_mu)}
% where eta is a small positive number to avoid division by zero.
% In case real(omega) == 0, eta is set to 0.0.
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

vc_data = isdftest.get(id_vc);
if isempty(vc_data.nrange1) || isempty(vc_data.nrange2)
  error('gen_Kq:nrange', ...
    'Missing cached nrange in ISDF id=%d. Build it first via isdftest.set_nrange(id, config.SYSTEM).', ...
    int32(id_vc));
end
nrangev = double(vc_data.nrange1);
nrangec = double(vc_data.nrange2);
Nisdf_vc = vc_data.nisdf;

Nisdf_o = Nisdf_vc;
R_sampling = vc_data.bundle_struct.R_grid_bundle(vc_data.bundle_struct.sampling2bundle, :);
fft_data = FFT.manager('get');
fftgrid = fft_data.fftgrid;
symm_data = symmetry.get();
inv_rot_index = symm_data.inv_rot_index;


iqrot = 1;
chiq_ISDF = zeros(Nisdf_o, Nisdf_o);
for ikbz = 1:nkbz
  ikibz = k_data.bz2ibz(ikbz);
  ikrot = k_data.bz2rot(ikbz);
  ik_qbz = r_lat_data.qindx_X(iqibz, ikbz, 1);
  iGo = r_lat_data.qindx_X(iqibz, ikbz, 2);
  ik_qibz = k_data.bz2ibz(ik_qbz);
  ik_qrot = k_data.bz2rot(ik_qbz);
  Go = double(r_lat_data.Ggrid_RLU(iGo, :));
  has_phase_shift = norm(Go) >= 1e-6;
  if has_phase_shift
    inviqrot = inv_rot_index(iqrot);
    invSqGo = double(Go * symm_data.rot_mtrx_RLU_G(:, :, inviqrot));
    phase_shift_coarse = exp(-2 * pi * 1i * (R_sampling ./ double(fftgrid)) * invSqGo');
  end
  nrangec_row = reshape(nrangec, 1, []);
  f_c = system_data.f(nrangec_row, ik_qibz, ispin);
  e_c = ev(nrangec_row, ik_qibz, ispin);
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
    is = [iv, ikibz, ikrot, ispin];
    u_xalpha_is = isdftest.get_u_xalpha(id_vc, is, iqrot);
    id_valid = find(valid);
    for it = 1:numel(id_valid)
      jc = nrangec_row(id_valid(it));
      os = [jc, ik_qibz, ik_qrot, ispin];
      u_xalpha_os = isdftest.get_u_xalpha(id_vc, os, iqrot);
      rho_col = conj(u_xalpha_is) .* u_xalpha_os;
      if has_phase_shift
        rho_col = rho_col .* phase_shift_coarse;
      end
      rho_blk(:, it) = rho_col;
    end
    rho_weighted = rho_blk .* reshape(coeff, 1, []);
    chiq_ISDF = chiq_ISDF + (rho_weighted * rho_blk');
  end
end % ikbz
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

% Then K
% K_q(\mu, \nu) = \Lambda_q^{-1}(\mu, \nu) - \tildeV_q(\mu, \nu)
Kq_ISDF = inv(chiq_ISDF) + vc_data.tildeVq(:, :, iqibz);

end