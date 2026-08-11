% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11

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

  isdf.free();
  report_dir = filename_map().isdf_report_dir;
  if exist(report_dir, 'dir') ~= 7
    mkdir(report_dir);
  end
  isdf.report.run_summary('clear');
  isdf.report.cond('init');

  if do_vc
    id_vc_new = run_isdf_type(config, cfg, 'vc');
    if cfg.validate_hf
      isdf.validation.isdf_validation(id_vc_new);
    end
  end

  if do_vn
    id_vn_new = run_isdf_type(config, cfg, 'vn');
    if cfg.validate_hf
      isdf.validation.isdf_validation(id_vn_new);
    end
  end

  if do_nn
    id_nn_new = run_isdf_type(config, cfg, 'nn');
    if cfg.validate_hf
      isdf.validation.isdf_validation(id_nn_new);
    end
  end

  isdf.report.run_summary('write');
end

function idnew = run_isdf_type(config, cfg, isdf_type)
  switch lower(isdf_type)
    case 'vc'
      id = isdf.isdf_add("vc");
      isdf_ratio = cfg.isdf_ratio_type1;
    case 'vn'
      id = isdf.isdf_add("vn");
      isdf_ratio = cfg.isdf_ratio_type2;
    case 'nn'
      id = isdf.isdf_add("nn");
      isdf_ratio = cfg.isdf_ratio_type3;
    otherwise
      error('isdf:driver:type', 'Unknown ISDF type ''%s''.', isdf_type);
  end

  isdf.set_nrange(id, config.SYSTEM);
  isdf_data = isdf.get(id);
  Nnrange1 = double(isdf_data.Nnrange1);
  Nnrange2 = double(isdf_data.Nnrange2);
  nmu_target = isdf_ratio * sqrt(Nnrange1 * Nnrange2);
  isdf_data.nisdf = int32(max(1, ceil(nmu_target)));
  isdf.save2mod(isdf_data, id);

  idx_mu = isdf.coeff.gen_coeff(cfg, id);
  if ~isempty(idx_mu)
    isdf.rsymm.init_from_indices(id, idx_mu);
  end
  isdf_data = isdf.get(id);
  isdf.report.run_summary('seed', struct( ...
    'desc', char(string(isdf_data.desc)), ...
    'seed_id', int32(id), ...
    'seed_nmu', int32(isdf_data.nisdf), ...
    'method', char(string(isdf_data.interp_scheme))));
  % isdf.coeff.print_coarse_grid_report(id);
  idnew = isdf.adaptive.launcher(id, cfg);
  local_dispatch_gen_tildeVq(idnew, config, cfg);
end

function local_dispatch_gen_tildeVq(id, config, cfg_isdf)
% Route to Gamma-optimized builder when frequency_dependence == -2.
  % if nargin >= 2 && isstruct(config) && isfield(config, 'FREQUENCY') ...
  %     && isfield(config.FREQUENCY, 'frequency_dependence') ...
  %     && config.FREQUENCY.frequency_dependence == -2
  %   fprintf('[ISDF] frequency_dependence=-2: gen_tildeVq_Gamma(id=%d)\n', int32(id));
  %   isdf.gen_tildeVq_Gamma(id, config);
  %   ratio = 0.75;
  %   if nargin >= 3 && isstruct(cfg_isdf) && isfield(cfg_isdf, 'inv_ratio')
  %     ratio = cfg_isdf.inv_ratio;
  %   end
  %   isdf.attach_gamma_trunc_factors(id, ratio);
  %   return;
  % end
  isdf.gen_tildeVq(id, cfg_isdf);
end
