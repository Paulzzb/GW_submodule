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
  isdf.report.cond('init');

  if do_vc
    id_vc_new = local_run_isdf_type(config, cfg, 'vc');
    if cfg.validate_hf
      isdf.validation.isdf_validation(id_vc_new);
    end
  end

  if do_vn
    id_vn_new = local_run_isdf_type(config, cfg, 'vn');
    if cfg.validate_hf
      isdf.validation.isdf_validation(id_vn_new);
    end
  end

  if do_nn
    id_nn_new = local_run_isdf_type(config, cfg, 'nn');
    if cfg.validate_hf
      isdf.validation.isdf_validation(id_nn_new);
    end
  end

end

function idnew = local_run_isdf_type(config, cfg, isdf_type)
  if local_use_sc_isdf(config)
    idnew = local_sc_isdf_pipeline(config, cfg, isdf_type);
    return;
  end

  switch lower(isdf_type)
    case 'vc'
      id = isdf.isdf_add("vc");
    case 'vn'
      id = isdf.isdf_add("vn");
    case 'nn'
      id = isdf.isdf_add("nn");
    otherwise
      error('isdf:driver:type', 'Unknown ISDF type ''%s''.', isdf_type);
  end

  isdf.set_nrange(id, config.SYSTEM);
  isdf_data = isdf.get(id);
  Nnrange1 = double(isdf_data.Nnrange1);
  Nnrange2 = double(isdf_data.Nnrange2);
  ratio_field = local_ratio_field(isdf_type);
  nmu_target = cfg.(ratio_field) * sqrt(Nnrange1 * Nnrange2);
  isdf_data.nisdf = int32(max(1, ceil(nmu_target)));
  isdf.save2mod(isdf_data, id);

  isdf.coeff.gen_coeff(cfg, id);
  isdf.coeff.print_coarse_grid_report(id);
  idnew = local_adaptiveisdf(id, cfg, isdf_type);
  local_dispatch_gen_tildeVq(idnew, config, cfg);
end

function tf = local_use_sc_isdf(config)
  tf = false;
  if ~isfield(config, 'SUPERCELL')
    return;
  end
  sc = config.SUPERCELL;
  if isfield(sc, 'use_sc_isdf') && logical(sc.use_sc_isdf)
    tf = true;
  end
end

function idnew = local_sc_isdf_pipeline(config, cfg, isdf_type)
  sc = config.SUPERCELL;
  ratio = int32([sc.k1, sc.k2, sc.k3]);
  source_dir = char(string(sc.isdf_source_dir));
  ck = local_sc_isdf_checkpoint(source_dir, isdf_type);
  fprintf('\n[ISDF SC_ISDF] type=%s source=%s ratio=[%d %d %d]\n', ...
    isdf_type, ck, ratio(1), ratio(2), ratio(3));
  ratio_field = local_ratio_field(isdf_type);
  target_ratio = cfg.(ratio_field);
  id_sc = isdf.SC_ISDF(ck, ratio, 'TargetIsdfRatio', target_ratio, ...
    'SystemCfg', config.SYSTEM);

  if local_sc_adaptive_enabled(sc)
    fprintf('[ISDF SC_ISDF] sc_adaptive=true: refine replicated grid with adaptiveisdf\n');
    isdf.set_nrange(id_sc, config.SYSTEM);
    isdf.SC_ISDF_prepare_adaptive_seed(id_sc);
    idnew = local_adaptiveisdf_sc(id_sc, cfg, isdf_type);
    isdf.get_u_xalpha('reset');
    local_dispatch_gen_tildeVq(idnew, config, cfg);
    return;
  end

  idnew = id_sc;
  isdf.coeff.gen_coeff_from_fine_grid(idnew);
  isdf.rsymm.gen_bundle(idnew);
  isdf.get_u_xalpha('reset');
  local_dispatch_gen_tildeVq(idnew, config, cfg);
end

function tf = local_sc_adaptive_enabled(sc)
  tf = false;
  if isfield(sc, 'sc_adaptive') && logical(sc.sc_adaptive)
    tf = true;
  end
end

function fpath = local_sc_isdf_checkpoint(source_dir, isdf_type)
  isdf_type = lower(strtrim(char(string(isdf_type))));
  patt = fullfile(source_dir, sprintf('isdf_adaptive_checkpoint_%s_id*.mat', isdf_type));
  d = dir(patt);
  if isempty(d)
    error('isdf:driver:SC_ISDF:MissingCheckpoint', ...
      'No checkpoint matching %s under %s.', patt, source_dir);
  end
  [~, ord] = sort([d.datenum], 'descend');
  fpath = fullfile(d(ord(1)).folder, d(ord(1)).name);
end

function field = local_ratio_field(isdf_type)
  switch lower(isdf_type)
    case 'vc'
      field = 'isdf_ratio_type1';
    case 'vn'
      field = 'isdf_ratio_type2';
    case 'nn'
      field = 'isdf_ratio_type3';
    otherwise
      field = 'isdf_ratio';
  end
end

function idnew = local_adaptiveisdf(id, cfg, isdf_type)
  cutoff = 1e-6;
  isdf_type = lower(strtrim(char(string(isdf_type))));
  isdf_data = isdf.get(id);
  desc = char(string(isdf_data.desc));

  switch isdf_type
    case 'vc'
      thr = local_adaptive_threshold(cfg, 'adaptive_threshold_type1');
      backend = 'adaptive_single';
      arithmetic = 'single isdf.prod, double Schur algebra';
      if thr < cutoff
        warning('isdf:driver:vcSinglePrecision', ...
          ['vc adaptive: threshold=%.4e < %.1e; ', ...
           'single-precision adaptive may be insufficient.'], thr, cutoff);
      end
      local_print_adaptive_dispatch(id, desc, isdf_type, backend, arithmetic, thr, cutoff);
      idnew = isdf.adaptive_single.adaptiveisdf(id, cfg);

    case 'vn'
      thr = local_adaptive_threshold(cfg, 'adaptive_threshold_type2');
      if thr > cutoff
        backend = 'adaptive_single';
        arithmetic = 'single isdf.prod, double Schur algebra';
      else
        backend = 'adaptive_double';
        arithmetic = 'double';
      end
      local_print_adaptive_dispatch(id, desc, isdf_type, backend, arithmetic, thr, cutoff);
      if thr > cutoff
        idnew = isdf.adaptive_single.adaptiveisdf(id, cfg);
      else
        idnew = isdf.adaptive_double.adaptiveisdf(id, cfg);
      end

    case 'nn'
      thr = local_adaptive_threshold(cfg, 'adaptive_threshold_type3');
      if thr > cutoff
        backend = 'adaptive_single';
        arithmetic = 'single isdf.prod, double Schur algebra';
      else
        backend = 'adaptive_double';
        arithmetic = 'double';
      end
      local_print_adaptive_dispatch(id, desc, isdf_type, backend, arithmetic, thr, cutoff);
      if thr > cutoff
        idnew = isdf.adaptive_single.adaptiveisdf(id, cfg);
      else
        idnew = isdf.adaptive_double.adaptiveisdf(id, cfg);
      end

    otherwise
      error('isdf:driver:adaptiveType', 'Unknown ISDF type ''%s''.', isdf_type);
  end
end

function local_print_adaptive_dispatch(id, desc, isdf_type, backend, arithmetic, thr, cutoff)
  fprintf(['\n[ISDF adaptive] id=%d desc=%s type=%s backend=%s\n', ...
    '  arithmetic=%s\n', ...
    '  threshold=%.8e (routing cutoff=%.1e)\n'], ...
    int32(id), desc, isdf_type, backend, arithmetic, thr, cutoff);
end

function idnew = local_adaptiveisdf_sc(id, cfg, isdf_type)
% SC adaptive refinement: always use double backend for numerical stability.
  isdf_type = lower(strtrim(char(string(isdf_type))));
  isdf_data = isdf.get(id);
  desc = char(string(isdf_data.desc));
  thr = local_adaptive_threshold(cfg, local_adaptive_threshold_field(isdf_type));
  local_print_adaptive_dispatch(id, desc, isdf_type, 'adaptive_double', ...
    'double (SC sc_adaptive)', thr, 1e-6);
  idnew = isdf.adaptive_double.adaptiveisdf(id, cfg);
end

function field = local_adaptive_threshold_field(isdf_type)
  switch lower(isdf_type)
    case 'vc'
      field = 'adaptive_threshold_type1';
    case 'vn'
      field = 'adaptive_threshold_type2';
    case 'nn'
      field = 'adaptive_threshold_type3';
    otherwise
      field = 'adaptive_threshold_type2';
  end
end

function thr = local_adaptive_threshold(cfg, field_name)
  thr = 2e-4;
  if isfield(cfg, field_name) && ~isempty(cfg.(field_name))
    v = double(cfg.(field_name));
    if isfinite(v) && v > 0
      thr = v;
    end
  end
end

function local_dispatch_gen_tildeVq(id, config, cfg_isdf)
% Route to Gamma-optimized builder when frequency_dependence == -2.
  if nargin >= 2 && isstruct(config) && isfield(config, 'FREQUENCY') ...
      && isfield(config.FREQUENCY, 'frequency_dependence') ...
      && config.FREQUENCY.frequency_dependence == -2
    fprintf('[ISDF] frequency_dependence=-2: testfunc(config, id=%d)\n', int32(id));
    isdf.testfunc(config, id);
    ratio = 0.75;
    if nargin >= 3 && isstruct(cfg_isdf) && isfield(cfg_isdf, 'inv_ratio')
      ratio = cfg_isdf.inv_ratio;
    end
    isdf.attach_gamma_trunc_factors(id, ratio);
    return;
  end
  isdf.gen_tildeVq(id, cfg_isdf);
end
