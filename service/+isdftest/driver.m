% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

function driver(~, config)

  cfg = config.ISDF;
  do_vc = true;
  do_vn = true;
  do_nn = true;
  if isfield(cfg, 'compute_vc')
    do_vc = logical(cfg.compute_vc);
  end
  if isfield(cfg, 'compute_vn')
    do_vn = logical(cfg.compute_vn);
  end
  if isfield(cfg, 'compute_nn')
    do_nn = logical(cfg.compute_nn);
  end

  % Keep in mind that,
  %     C_{n1k1n2k2, ga} = \conj{coeff_seperated(n1k1, ga)} * coeff_seperated(n2k2, ga)
  % coeff_seperated: size nisdf * nb * nkibz * nspin.

  isdftest.free();
  isdftest.numerical_cond_report('init');
  if do_vc
    id_vc = isdftest.isdftest_add("vc");
    isdftest.set_nrange(id_vc, config.SYSTEM);
    vc_data = isdftest.get(id_vc);
    Nnrange1 = double(vc_data.Nnrange1);
    Nnrange2 = double(vc_data.Nnrange2);
    nmu_target = cfg.isdf_ratio_type1 * sqrt(Nnrange1 * Nnrange2);
    vc_data.nisdf = int32(max(1, ceil(nmu_target)));
    isdftest.save2mod(vc_data, id_vc);
    %
    % Coarse-grid ISDF pipeline (ignore config.isdftest.exxmethod stubs such as kmeans).
    isdftest.coeff.gen_coeff(cfg, id_vc);
    isdftest.coeff.print_coarse_grid_report(id_vc);
    %
    id_vc_new = isdftest.adaptive.adaptiveisdf(id_vc, cfg);
    %
    isdftest.gen_tildeVq(id_vc_new, cfg);
    %
    isdftest.validation.isdf_validation(id_vc_new);
  end

  if do_vn
    id_vn = isdftest.isdftest_add("vn");
    isdftest.set_nrange(id_vn, config.SYSTEM);
    vn_data = isdftest.get(id_vn);
    Nnrange1 = double(vn_data.Nnrange1);
    Nnrange2 = double(vn_data.Nnrange2);
    nmu_target = cfg.isdf_ratio_type2 * sqrt(Nnrange1 * Nnrange2);
    vn_data.nisdf = int32(max(1, ceil(nmu_target)));
    isdftest.save2mod(vn_data, id_vn);
    %
    % Coarse-grid ISDF pipeline (ignore config.isdftest.exxmethod stubs such as kmeans).
    isdftest.coeff.gen_coeff(cfg, id_vn);
    isdftest.coeff.print_coarse_grid_report(id_vn);
    %
    id_vn_new = isdftest.adaptive.adaptiveisdf(id_vn, cfg);
    %
    isdftest.gen_tildeVq(id_vn_new, cfg);
    %
    isdftest.validation.isdf_validation(id_vn_new);
  end

  if do_nn
    id_nn = isdftest.isdftest_add("nn");
    isdftest.set_nrange(id_nn, config.SYSTEM);
    nn_data = isdftest.get(id_nn);
    Nnrange1 = double(nn_data.Nnrange1);
    Nnrange2 = double(nn_data.Nnrange2);
    nmu_target = cfg.isdf_ratio_type3 * sqrt(Nnrange1 * Nnrange2);
    nn_data.nisdf = int32(max(1, ceil(nmu_target)));
    % nn_data.assigned = true;
    isdftest.save2mod(nn_data, id_nn);
    %
    isdftest.coeff.gen_coeff(cfg, id_nn);
    isdftest.coeff.print_coarse_grid_report(id_nn);
    %
    id_nn_new = isdftest.adaptive.adaptiveisdf(id_nn, cfg);
    %
    isdftest.gen_tildeVq(id_nn_new, cfg);
    %
    isdftest.validation.isdf_validation(id_nn_new);
  end



end
