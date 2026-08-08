% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ
function Kq_ISDF = gen_Kq_Gamma(id_vc, iqibz, omega, broadening, cauchy)
%GEN_KQ_GAMMA  Build K_q at Gamma (single-k BZ).
%
%   Kq_ISDF = isdf.gen_Kq_Gamma(id_vc, iqibz)
%   Kq_ISDF = isdf.gen_Kq_Gamma(id_vc, iqibz, omega)
%   Kq_ISDF = isdf.gen_Kq_Gamma(id_vc, iqibz, omega, broadening)
%   Kq_ISDF = isdf.gen_Kq_Gamma(id_vc, iqibz, omega, broadening, cauchy)
%
%   id_vc       — ISDF slot for polarizability (typically desc 'vc')
%   iqibz       — q IBZ index
%   omega       — frequency in Ry (same as system.Eo); default 0
%   broadening  — eta in Ry (caller converts FREQUENCY.broadening eV → Ry);
%                 omit / [] → eta = 0
%   cauchy      — [] | struct from isdf.cauchy_opts (isCauchy/froErr/MaxIter)
%
% Chi path:
%   eta ~= 0  → always direct iv-loop (real-axis)
%   eta == 0  → if cauchy.isCauchy, use isdf.Cauchy.COmegaCstar (static Ω);
%               otherwise direct loop. Imag-axis fullfreq passes broadening=[]
%               so Cauchy is eligible when enabled.
%
%   K_q = chi_q^{-1} + V~_q   (after SVD truncation factors on chi)

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

if nargin < 3 || isempty(omega)
  omega = complex(0);
end
if nargin < 4 || isempty(broadening)
  eta = 0;
else
  eta = double(broadening);
end
if nargin < 5 || isempty(cauchy)
  cauchy = struct('isCauchy', false);
end

nspin = 1; ispin = 1;

system_data = system.get();
k_data = lattice.manager('k', 'get');
nkbz = k_data.nbz;
if nkbz ~= 1
  output.err('Only single-k point sampling is supported.');
end
ev = system_data.Eo;

vc_data = isdf.get(id_vc);
if isempty(vc_data.nrange1) || isempty(vc_data.nrange2)
  output.err( ...
    'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', ...
    int32(id_vc));
end
nrangev = double(vc_data.nrange1);
nrangec = double(vc_data.nrange2);
Nisdf_vc = vc_data.nisdf;
s2b = vc_data.bundle_struct.sampling2bundle;
Nisdf_o = Nisdf_vc;

iqrot = 1;
ikibz = 1;
nrangec_row = reshape(nrangec, 1, []);

use_cauchy = (abs(eta) < 1e-14) && isfield(cauchy, 'isCauchy') ...
  && logical(cauchy.isCauchy);

if use_cauchy
  % COmegaCstar is static (Ω_ij = e_v - e_c). Nonzero omega → direct path.
  if abs(omega) > 1e-12
    output.msg('v1s', ...
      'Cauchy requested but omega~=0; using direct chi sum.');
    use_cauchy = false;
  end
end

if use_cauchy
  Phi = vc_data.bundle_struct.WF_bundle(s2b, nrangev, 1, ispin);
  Psi = vc_data.bundle_struct.WF_bundle(s2b, nrangec_row, 1, ispin);
  evOcc = double(ev(nrangev, ikibz, ispin));
  evUnocc = double(ev(nrangec_row, ikibz, ispin));
  optC = struct('froErr', 1e-6, 'MaxIter', 10);
  if isfield(cauchy, 'froErr') && ~isempty(cauchy.froErr)
    optC.froErr = double(cauchy.froErr);
  end
  if isfield(cauchy, 'MaxIter') && ~isempty(cauchy.MaxIter)
    optC.MaxIter = double(cauchy.MaxIter);
  end
  [chi_raw, ~, ~] = isdf.Cauchy.COmegaCstar(Phi, Psi, evOcc(:), evUnocc(:), optC);
  % Static limit of direct coeff = occ/den1-occ/den2 is -2/(Ev-Ec).
  % COmegaCstar builds 1/(Ev-Ec); negate and *2 match that (then spin *2).
  chiq_ISDF = -chi_raw;
  chiq_ISDF = 2.0 * chiq_ISDF;
  if nspin == 1
    chiq_ISDF = 2 * chiq_ISDF;
  end
else
  chiq_ISDF = zeros(Nisdf_o, Nisdf_o);
  f_c = system_data.f(nrangec_row, 1, ispin);
  e_c = ev(nrangec_row, 1, ispin);
  for iv = nrangev
    f_ik = system_data.f(iv, ikibz, ispin);
    e_ik = ev(iv, ikibz, ispin);
    occ = -f_c + f_ik;
    Delta = e_ik - e_c;
    if abs(eta) < 1e-14
      den1 = omega - Delta;
      den2 = omega + Delta;
    else
      den1 = omega - Delta - 1i * eta;
      den2 = omega + Delta + 1i * eta;
    end
    valid = (abs(occ) >= 1e-5) & (abs(den1) >= 1e-12) & (abs(den2) >= 1e-12);
    if ~any(valid)
      continue;
    end

    coeff = occ(valid) ./ den1(valid) - occ(valid) ./ den2(valid);
    is = [iv, 1, 1, ispin];
    u_xalpha_is = isdf.get_u_xalpha(id_vc, is, iqrot);
    u_xalpha_os = vc_data.bundle_struct.WF_bundle(s2b, nrangec_row, 1, ispin);
    rho_blk = conj(u_xalpha_is) .* u_xalpha_os(:, valid);
    rho_weighted = rho_blk .* reshape(coeff, 1, []);
    chiq_ISDF = chiq_ISDF + (rho_weighted * rho_blk');
  end
  % Factor 2: spin-unpolarized (nspin == 1); both poles already in coeff.
  if nspin == 1
    chiq_ISDF = 2 * chiq_ISDF;
  end
end

% Hermitianize when eta = 0 (COHSEX / imag-axis)
if abs(eta) < 1e-6
  chiq_ISDF = (chiq_ISDF + chiq_ISDF') / 2;
end

% Apply Lambda^{inv_ratio-1} V truncations
inv_ratio = vc_data.svd_ratio;
fac = vc_data.CCHq_trunc_factors{double(iqibz)};
V_trunc = fac.V_trunc;
Lambda_trunc = fac.Lambda_trunc;
N_keep = fac.N_keep;

chiq_ISDF = chiq_ISDF * V_trunc;
chiq_ISDF = V_trunc' * chiq_ISDF;
chiq_ISDF = Lambda_trunc.^(inv_ratio-1) .* (chiq_ISDF) .* Lambda_trunc.^(inv_ratio-1).';
if abs(eta) < 1e-6
  chiq_ISDF = 0.5 * chiq_ISDF' + 0.5 * chiq_ISDF;
end

% K_q = chi_q^{-1} + V~_q
Kq_ISDF = inv(chiq_ISDF) + vc_data.tildeVq(1:N_keep, 1:N_keep, iqibz);

end
