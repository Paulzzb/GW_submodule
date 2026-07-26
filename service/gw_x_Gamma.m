% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23 ZZ

function Ex = gw_x_Gamma(config)

msg = sprintf('[Exchange] Start computing Sigma_x (exchange part) with service packages...\n');
QPlog(msg, 0);
tStart = tic;

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

% Basic setup from GWinfo/config
nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
nspin = 1;
warning('gw_x_k_packages: nspin and ispin are set to 1 for now.');

% Pull runtime package data from service managers
r_lat_m = lattice.manager('r_lat', 'get');
k = lattice.manager('k', 'get');
q = lattice.manager('q', 'get');
coul_data = coulomb.manager('get');
system_data = system.get();
% wf_data = wf.get();

nb = system_data.nb;

% Key arrays for Sigma_x loop
qindx_S = r_lat_m.qindx_S;
nkibz = k.nibz;
nqibz = q.nibz;
nqbz = q.nbz;
sstar = [k.bz2ibz, k.bz2rot];
vcoul_all = coul_data.vcoul;

if isempty(vcoul_all)
  error('gw_x_k_packages: vcoul is empty. Run coulomb.driver (and coulomb.vcoul) first.');
end

if size(vcoul_all, 2) < nqibz
  error('gw_x_k_packages: vcoul column size (%d) is smaller than nqibz (%d).', size(vcoul_all, 2), nqibz);
end

msg = sprintf('[Exchange] Using service-package Sigma_x calculation.\n');
QPlog(msg, 0);

Ex = zeros(nbmax - nbmin + 1, nkibz, nspin);
tStandard = tic;

isisdf = logical(config.ISDF.isisdf);
id_ex = [];
ex_data = [];
if isisdf
  [id_ex, ex_label] = local_resolve_ex_isdf_slot(config);
  QPlog(sprintf('[Exchange] ISDF Sigma_x: %s adaptive slot + get_rho_xalpha + tildeVq.', ex_label), 0);
  if isempty(id_ex)
    error('gw_x_k_packages:isdfExSlot', ...
      ['ISDF exchange requires %s slot (adaptive). ' ...
       'Enable config.ISDF.compute_%s or run isdf.driver first.'], ex_label, ex_label);
  end
  isdf.set_nrange(id_ex, config.SYSTEM);
  ex_data = isdf.get(id_ex);
  if isempty(ex_data.tildeVq) || size(ex_data.tildeVq, 3) < nqibz
    QPlog(sprintf('[Exchange] Building tildeVq for %s id=%d ...', ex_label, id_ex), 1);
    isdf.gen_tildeVq(id_ex, config.ISDF);
    ex_data = isdf.get(id_ex);
  end
  fprintf('[Exchange] ISDF %s id=%d, nisdf=%d, scheme=%s\n', ...
    ex_label, id_ex, ex_data.nisdf, char(ex_data.interp_scheme));
  s2b_ex = ex_data.bundle_struct.sampling2bundle;
  fac_ex = ex_data.CCHq_trunc_factors{1};
  s_ratio_ex = ex_data.svd_ratio;
  Nkeep_ex = fac_ex.N_keep;
else
  QPlog('[Exchange] Dense G-space Sigma_x (SCATTER_Bamp).', 0);
end

for ik = 1:nkibz
% for ik = 1:nkibz
  for ib = nbmin:nbmax
    ib_out = ib - nbmin + 1;
    for ispin = 1:nspin
      for iq = 1:nqbz
        iqibz = sstar(iq, 1);
        iqs = sstar(iq, 2);

        ikp_bz = qindx_S(ik, iq, 1);
        is = qindx_S(ik, iq, 2);
        ikp_ibz = sstar(ikp_bz, 1);
        ikp_rot = sstar(ikp_bz, 2);

        vcoul_q = vcoul_all(:, iqibz);
        if iqibz == 1
          vcoul_q(1) = coul_data.vcoul0; % Set G=0 component to zero for exchange term
        end

        if isisdf
          u_ib_xalpha = ex_data.bundle_struct.WF_bundle(s2b_ex, ib, 1, ispin);
          u_ob_xalpha = zeros(ex_data.nisdf, nb);
          max_occ = 0;
          for ob = 1:nb 
            occ = system_data.f(ob, ikp_ibz, ispin);
            if occ < 1e-6
              continue;
            end
            max_occ = ob;
            % Same ISDF contraction as isdf_validate_HF / gw_cohsex (adaptive + tildeVq).
            u_ob_xalpha(:, ob) = ex_data.bundle_struct.WF_bundle(s2b_ex, ob, 1, ispin);
          end

          rho_left_ex = conj(u_ib_xalpha) .* u_ob_xalpha(:, 1:max_occ);
          rho_left_ex = ...
            diag(fac_ex.Lambda_trunc.^(s_ratio_ex-1)) * (fac_ex.V_trunc' * rho_left_ex);
          tmp1 = ex_data.tildeVq(1:Nkeep_ex, 1:Nkeep_ex, 1) * rho_left_ex;
          ex_t = sum(conj(rho_left_ex) .* tmp1, 1);
          ex_t = real(ex_t);
          Ex_t = sum(ex_t);

          Ex(ib_out, ik, ispin) = Ex(ib_out, ik, ispin) - Ex_t;
        else
          for ob = 1:nb 
            occ = system_data.f(ob, ikp_ibz, ispin);
            if occ < 1e-6
              continue;
            end

            param = [];
            param.is = [ib, ik, 1, ispin];
            param.os = [ob, ikp_ibz, ikp_rot, ispin];
            param.qs = [is, iqibz, iqs];

            if isisdf
              % Same ISDF contraction as isdf_validate_HF / gw_cohsex (adaptive + tildeVq).
              rho_mu = isdf.get_rho_xalpha(id_ex, param);
              tildeVq_q = ex_data.tildeVq(:, :, iqibz);
              Ex_t = real(rho_mu' * tildeVq_q * rho_mu);
              ngrho_left = SCATTER_Bamp(param);
              Ex_t_dir = sum(vcoul_q .* abs(ngrho_left).^2);
              fprintf('param: is=[%d %d %d %d] os=[%d %d %d %d] qs=[%d %d %d]\n', ...
                param.is(1), param.is(2), param.is(3), param.is(4), ...
                param.os(1), param.os(2), param.os(3), param.os(4), ...
                param.qs(1), param.qs(2), param.qs(3));
              fprintf('  Ex_t=%.6f  Ex_t_dir=%.6f  diff=%.6f\n', Ex_t, Ex_t_dir, Ex_t - Ex_t_dir);
            else
              ngrho_left = SCATTER_Bamp(param);
              Ex_t = sum(vcoul_q .* abs(ngrho_left).^2);
            end
            Ex(ib_out, ik, ispin) = Ex(ib_out, ik, ispin) - Ex_t;
          end
        end
      end
    end
  end
end

msg = sprintf('[Exchange] Standard loop completed in %.2f seconds.\n', toc(tStandard));
QPlog(msg, 1);

msg = sprintf('[Exchange] Finished. Total time: %.2f seconds.\n', toc(tStart));
QPlog(msg, 0);

end % EOF

function [id_ex, ex_label] = local_resolve_ex_isdf_slot(config)
  ex_label = "vn";
  if isfield(config, 'COHSEX') && isfield(config.COHSEX, 'ex_use_which_isdf') ...
      && ~isempty(config.COHSEX.ex_use_which_isdf)
    ex_label = lower(string(config.COHSEX.ex_use_which_isdf));
  end
  if ~ismember(ex_label, ["vn", "nn"])
    error('gw_x_k_packages:ExUseWhichIsdf', ...
      'COHSEX.ex_use_which_isdf must be ''vn'' or ''nn'', got ''%s''.', ex_label);
  end
  [~, id_vn, id_nn] = isdf.cohsex_resolve_ids(config);
  if ex_label == "vn"
    id_ex = id_vn;
  else
    id_ex = id_nn;
  end
end
