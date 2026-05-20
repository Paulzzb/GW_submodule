% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23 ZZ

function Ex = gw_x_k_packages(config)

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
ispin = 1;
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

          ngrho_left = SCATTER_Bamp(param);

          Ex_t = sum(vcoul_q .* abs(ngrho_left).^2);
          Ex(ib_out, ik, ispin) = Ex(ib_out, ik, ispin) - Ex_t;
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
