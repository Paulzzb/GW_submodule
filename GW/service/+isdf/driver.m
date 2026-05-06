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

  isdf.free();
  if do_vc
    id_vc = isdf.isdf_add("vc");
    [nrange1, nrange2] = isdf.isdf_get_nrange(id_vc);
    vc_data = isdf.get(id_vc);
    nmu_target = cfg.isdf_ratio_type1 * sqrt(length(nrange1) * length(nrange2));
    vc_data.nisdf = int32(max(1, ceil(nmu_target)));
    isdf.save2mod(vc_data, id_vc);
    %
    % Coarse-grid ISDF pipeline (ignore config.ISDF.exxmethod stubs such as kmeans).
    isdf.coeff.gen_coeff(cfg, id_vc);
    isdf.coeff.print_coarse_grid_report(id_vc);
    %
    id_vc_new = isdf.adaptive.adaptiveisdf(id_vc, cfg);
    %
    isdf.gen_tildeVq(id_vc_new);
    %
    isdf.validation.isdf_validation(id_vc_new);
  end

  if do_vn
    id_vn = isdf.isdf_add("vn");
    [nrange1, nrange2] = isdf.isdf_get_nrange(id_vn);
    vn_data = isdf.get(id_vn);
    nmu_target = cfg.isdf_ratio_type2 * sqrt(length(nrange1) * length(nrange2));
    vn_data.nisdf = int32(max(1, ceil(nmu_target)));
    isdf.save2mod(vn_data, id_vn);
    %
    % Coarse-grid ISDF pipeline (ignore config.ISDF.exxmethod stubs such as kmeans).
    isdf.coeff.gen_coeff(cfg, id_vn);
    isdf.coeff.print_coarse_grid_report(id_vn);
    %
    id_vn_new = isdf.adaptive.adaptiveisdf(id_vn, cfg);
    %
    isdf.gen_tildeVq(id_vn_new);
    %
    isdf.validation.isdf_validation(id_vn_new);
  end

  if do_nn
    id_nn = isdf.isdf_add("nn");
    nn_data = isdf.get(id_nn);
    [nrange1, nrange2] = isdf.isdf_get_nrange(id_nn);
    nmu_target = cfg.isdf_ratio_type3 * sqrt(length(nrange1) * length(nrange2));
    nn_data.nisdf = int32(max(1, ceil(nmu_target)));
    % nn_data.assigned = true;
    isdf.save2mod(nn_data, id_nn);
    %
    isdf.coeff.gen_coeff(cfg, id_nn);
    isdf.coeff.print_coarse_grid_report(id_nn);
    %
    id_nn_new = isdf.adaptive.adaptiveisdf(id_nn, cfg);
    %
    isdf.gen_tildeVq(id_nn_new);
    %
    isdf.validation.isdf_validation(id_nn_new);
  end



end
