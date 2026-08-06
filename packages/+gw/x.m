% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function Ex = x(config)

msg = sprintf('[Exchange] Start computing Sigma_x (exchange part) with service packages...\n');
output.msg('v0s', '%s', msg);
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
msg = sprintf('nspin and ispin are set to 1 for now.\n');
output.warn('%s', msg);

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
  output.err('vcoul is empty. Run coulomb.driver (and coulomb.vcoul) first.');
end

if size(vcoul_all, 2) < nqibz
  output.err('vcoul column size (%d) is smaller than nqibz (%d).', size(vcoul_all, 2), nqibz);
end

msg = sprintf('[Exchange] Using service-package Sigma_x calculation.\n');
output.msg('v0s', '%s', msg);

Ex = zeros(nbmax - nbmin + 1, nkibz, nspin);
tStandard = tic;

isisdf = logical(config.ISDF.isisdf);
id_ex = [];
ex_data = [];
if isisdf
  [id_ex, ex_label] = local_resolve_ex_isdf_slot(config);
  msg = sprintf('[Exchange] ISDF Sigma_x: %s adaptive slot + get_rho_xalpha + tildeVq.\n', ex_label);
  output.msg('v0s', '%s', msg);
  if isempty(id_ex)
    output.err( ...
      ['ISDF exchange requires %s slot (adaptive). ' ...
       'Enable config.ISDF.compute_%s or run isdf.driver first.'], ex_label, ex_label);
  end
  isdf.set_nrange(id_ex, config.SYSTEM);
  ex_data = isdf.get(id_ex);
  if isempty(ex_data.tildeVq) || size(ex_data.tildeVq, 3) < nqibz
    msg = sprintf('[Exchange] Building tildeVq for %s id=%d ...\n', ex_label, id_ex);
    output.msg('v1s', '%s', msg);
    isdf.gen_tildeVq(id_ex, config.ISDF);
    ex_data = isdf.get(id_ex);
  end
  msg = sprintf('[Exchange] ISDF %s id=%d, nisdf=%d, scheme=%s\n', ...
    ex_label, id_ex, ex_data.nisdf, char(ex_data.interp_scheme));
  output.msg('v0s', '%s', msg);
else
  msg = sprintf('[Exchange] Dense G-space Sigma_x (SCATTER_Bamp).\n');
  output.msg('v0s', '%s', msg);
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

msg = sprintf('[Exchange] Standard loop completed in %.2f seconds.\n', toc(tStandard));
output.msg('v1s', '%s', msg);

msg = sprintf('[Exchange] Finished. Total time: %.2f seconds.\n', toc(tStart));
output.msg('v0s', '%s', msg);

end % EOF

function [id_ex, ex_label] = local_resolve_ex_isdf_slot(config)
  ex_label = "vn";
  if isfield(config, 'COHSEX') && isfield(config.COHSEX, 'ex_use_which_isdf') ...
      && ~isempty(config.COHSEX.ex_use_which_isdf)
    ex_label = lower(string(config.COHSEX.ex_use_which_isdf));
  end
  if ~ismember(ex_label, ["vn", "nn"])
    output.err( ...
      'COHSEX.ex_use_which_isdf must be ''vn'' or ''nn'', got ''%s''.', ex_label);
  end
  [~, id_vn, id_nn] = isdf.cohsex_resolve_ids(config);
  if ex_label == "vn"
    id_ex = id_vn;
  else
    id_ex = id_nn;
  end
end
