% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

function driver(~, config)

  system_data = system.get();
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  nb = wf_data.nb;
  nkibz = k_data.nibz;

  cfg = config.ISDF;
  nocc_max = 0;
  nspin = system_data.nspin;
  for ispin = 1:nspin
    for ikibz = 1:nkibz
      f_ib = system_data.f(:, ikibz, ispin);
      idx_last = find(f_ib(:) > 1e-5, 1, 'last');
      if ~isempty(idx_last)
        nocc_max = max(nocc_max, idx_last);
      end
    end
  end

  % Keep in mind that,
  %     C_{n1k1n2k2, ga} = \conj{coeff_seperated(n1k1, ga)} * coeff_seperated(n2k2, ga)
  % coeff_seperated: size nisdf * nb * nkibz * nspin.
  
  isdf.free();
  % id_vn = isdf.isdf_add("vn");
  % vn_data = isdf.get(id_vn);
  % nmu_target = cfg.isdf_ratio_type2 * sqrt(double(nocc_max) * double(nb));
  % vn_data.nisdf = int32(max(1, ceil(nmu_target)));
  % isdf.save2mod(vn_data, id_vn);
  % %
  % % Coarse-grid ISDF pipeline (ignore config.ISDF.exxmethod stubs such as kmeans).
  % isdf.gen_coeff(cfg, id_vn);
  % isdf.print_coarse_grid_report(id_vn);
  % %
  % isdf.gen_tildeVq(id_vn);
  % %
  % isdf.isdf_validation(id_vn);
  


  id_nn = isdf.isdf_add("nn");
  nn_data = isdf.get(id_nn);
  nmu_target = cfg.isdf_ratio_type3 * double(nb);
  nn_data.nisdf = int32(max(1, ceil(nmu_target)));
  nn_data.assigned = true;
  isdf.save2mod(nn_data, id_nn);
  %
  isdf.gen_coeff(cfg, id_nn);
  isdf.print_coarse_grid_report(id_nn);
  %
  isdf.gen_tildeVq(id_nn);
  %
  isdf.isdf_validation(id_nn);



end
