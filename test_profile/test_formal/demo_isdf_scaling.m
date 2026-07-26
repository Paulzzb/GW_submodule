function report = demo_isdf_scaling(data, isdf_type)
%DEMO_ISDF_SCALING  Formal SC ISDF for vc + nn (no isdf.driver).
%
%   report = demo_isdf_scaling(data)
%
%   Always builds vc and nn.  Stores raw helperqG (before /CCHq), then forms
%   vcVvc and vcVnn via G-space Coulomb contraction followed by CCHq division:
%     M_ab = helperqG_a' * diag(vcoul) * helperqG_b
%     V_ab = (CCHq_a \ M_ab) / CCHq_b
%
%   Cache: UC_save_dir/demo_isdf_scaling_k<k1>_<k2>_<k3>.mat

  if nargin < 1 || ~isstruct(data)
    error('demo_isdf_scaling:data', 'demo_sc_scaling output struct is required.');
  end
  req = {'data_sc', 'config_sc', 'UC_save_dir', 'k1k2k3', 'data_uc'};
  for k = 1:numel(req)
    if ~isfield(data, req{k})
      error('demo_isdf_scaling:field', 'Missing field data.%s.', req{k});
    end
  end

  data_sc = data.data_sc;
  config = data.config_sc;
  UC_save_dir = char(string(data.UC_save_dir));
  ratio = int32(data.k1k2k3(:)).';
  uc_fftgrid = int32([data.data_uc.sys.n1, data.data_uc.sys.n2, data.data_uc.sys.n3]);

  if nargin >= 2 && ~isempty(isdf_type)
    types = {lower(strtrim(char(string(isdf_type))))};
  else
    types = {'vc', 'nn'};
  end

  batch_bytes = 512 * 1024^3;
  use_cache = true;
  if isfield(config, 'TESTFUNC') && isstruct(config.TESTFUNC)
    tf = config.TESTFUNC;
    if isfield(tf, 'mchq_batch_mib') && ~isempty(tf.mchq_batch_mib)
      batch_bytes = double(tf.mchq_batch_mib) * 1024^2;
    end
    if isfield(tf, 'use_cache') && ~isempty(tf.use_cache)
      use_cache = logical(tf.use_cache);
    end
  end

  cache_path = local_cache_path(UC_save_dir, types, ratio);
  if use_cache
    report = local_try_load_cache(cache_path, UC_save_dir, types, ratio, config, data_sc);
    if ~isempty(report)
      return
    end
  end

  local_bootstrap_service(data_sc, config);

  report = struct('types', string(types), 'ratio', double(ratio), ...
    'UC_save_dir', UC_save_dir, 'from_cache', false, 'cache_path', cache_path);

  for it = 1:numel(types)
    typ = types{it};
    ck = local_checkpoint(UC_save_dir, typ);
    fprintf('\n[demo_isdf_scaling] type=%s ck=%s ratio=[%d %d %d]\n', ...
      typ, ck, ratio(1), ratio(2), ratio(3));

    target_ratio = local_target_isdf_ratio(config, typ);
    id = isdf.SC_ISDF(ck, ratio, 'UcFftgrid', uc_fftgrid, ...
      'TargetIsdfRatio', target_ratio, 'SystemCfg', config.SYSTEM);
    isdf.set_nrange(id, config.SYSTEM);
    isdf.coeff.gen_coeff_from_fine_grid(id);
    isdf.rsymm.gen_bundle(id);
    isdf.get_u_xalpha('reset');

    case_out = local_build_helperqG(id, batch_bytes);
    case_out.isdf_type = string(typ);
    case_out.id = int32(id);
    case_out.ck_file = string(ck);
    report.(typ) = case_out;

    if strcmp(typ, 'vc')
      report.vcVvc = local_mu_coul_contract(case_out, case_out);
      report.vc.tildeVq = report.vcVvc;
      fprintf('[demo_isdf_scaling] vcVvc %dx%d\n', ...
        size(report.vcVvc, 1), size(report.vcVvc, 2));
    elseif strcmp(typ, 'nn')
      if ~isfield(report, 'vc')
        error('demo_isdf_scaling:order', 'Build vc before nn.');
      end
      report.vcVnn = local_mu_coul_contract(report.vc, case_out);
      fprintf('[demo_isdf_scaling] vcVnn %dx%d\n', ...
        size(report.vcVnn, 1), size(report.vcVnn, 2));
    end
  end

  if use_cache
    cache_meta = local_cache_meta(UC_save_dir, types, ratio, config, data_sc);
    save(cache_path, 'report', 'cache_meta', '-v7.3');
    fprintf('[demo_isdf_scaling] saved cache: %s (types: %s)\n', ...
      cache_path, strjoin(types, ', '));
  end
end

%% --- compute ---

function local_bootstrap_service(data_sc, config)
  test_profile_dir = fileparts(mfilename('fullpath'));
  gw_root = fileparts(fileparts(test_profile_dir));
  here = pwd;
  cleanup = onCleanup(@() cd(here)); %#ok<NASGU>

  cd(fullfile(gw_root, 'service'));
  addpath(genpath(pwd));
  rehash;
  local_ensure_mex();

  cd(gw_root);
  QPstartup;
  service_reset_persistent(struct('keep_pool', true));
  packages_reset_persistent();
  isdf.free();

  if ~isfield(config, 'FREQUENCY') || ~isstruct(config.FREQUENCY)
    config.FREQUENCY = default_param_values().FREQUENCY;
  end
  config.FREQUENCY.frequency_dependence = -2;
  config = set_default_param_value(config, data_sc);

  isdf.debug.init_from_config(config);
  parallel.driver(data_sc, config);
  system.driver(data_sc, config);
  symmetry.driver(data_sc, config);
  FFT.driver(data_sc, config);
  lattice.driver(data_sc, config);
  coulomb.driver(data_sc, config);
  wave_functions.driver(data_sc, config);
end

function case_out = local_build_helperqG(id, batch_bytes)
  wf = wave_functions.get();
  coul = coulomb.get();
  fft = FFT.get();
  lat = lattice.manager('d_lat', 'get');
  isdf_data = isdf.get(id);

  nrange1 = double(isdf_data.nrange1);
  nrange2 = double(isdf_data.nrange2);
  fft_sz = fft.fftgrid;
  g_col = fft.G_table(:, 1);
  DL_vol = lat.DL_vol;
  nc = double(wf.nc);
  ng = double(coul.coulomb_ng);
  Nmu = double(isdf_data.bundle_struct.N_sampling);
  n_batch = max(1, min(Nmu, floor(batch_bytes / (nc * 8))));

  fprintf(['[demo_isdf_scaling] id=%d desc=%s | nc=%.3g ng=%d Nmu=%.0f ', ...
    'nb1=%d nb2=%d | batch=%d\n'], int32(id), char(string(isdf_data.desc)), ...
    nc, ng, Nmu, numel(nrange1), numel(nrange2), n_batch);

  vcoul_q = coul.vcoul(:, 1);
  vcoul_q(1) = coul.vcoul0;

  if isfield(isdf_data.bundle_struct, 'fine_grid_lin') ...
      && ~isempty(isdf_data.bundle_struct.fine_grid_lin)
    fine_lin = int32(isdf_data.bundle_struct.fine_grid_lin(:));
  else
    fine_lin = isdf.SC_ISDF_r_sampling_to_lin(isdf_data);
  end

  wf_c1 = conj(wf.c(fine_lin, nrange1, 1, 1));
  wf_c2 = wf.c(fine_lin, nrange2, 1, 1);
  psixga = wf.c(fine_lin, :, 1, 1);
  CCHq = isdf.prod(wf_c1, wf_c1, wf_c2, wf_c2);
  CCHq = 0.5 * (CCHq + CCHq');

  wf1 = conj(wf.c(:, nrange1, 1, 1));
  wf2 = wf.c(:, nrange2, 1, 1);
  helperqG = zeros(ng, Nmu);

  for i0 = 1:n_batch:Nmu
    i1 = min(i0 + n_batch - 1, Nmu);
    idx = i0:i1;
    MCHq_R = isdf.prod(wf1, wf_c1(idx, :), wf2, wf_c2(idx, :));
    nbat = size(MCHq_R, 2);
    batch_G = zeros(ng, nbat);
    parfor ib = 1:nbat
      box = reshape(MCHq_R(:, ib), fft_sz);
      box = do_FFT(box, fft_sz, 1) * DL_vol;
      batch_G(:, ib) = box(g_col);
    end
    helperqG(:, idx) = batch_G;
  end

  case_out = struct('desc', string(isdf_data.desc), 'Nmu', Nmu, 'ng', ng, ...
    'nc', nc, 'n_batch', n_batch, 'helperqG', helperqG, 'CCHq', CCHq, ...
    'vcoul_q', vcoul_q, 'fine_grid_lin', fine_lin, 'psixga', psixga, ...
    'nrange1', nrange1, 'nrange2', nrange2);
end

function V = local_mu_coul_contract(rep_a, rep_b)
  D = rep_a.vcoul_q(:);
  M = rep_a.helperqG' * (D .* rep_b.helperqG);
  V = (rep_a.CCHq \ M) / rep_b.CCHq;
  if size(V, 1) == size(V, 2)
    V = 0.5 * (V + V');
  end
end

%% --- cache ---

function report = local_try_load_cache(cache_path, UC_save_dir, types, ratio, config, data_sc)
  report = [];
  if ~isfile(cache_path)
    return
  end
  S = load(cache_path, 'report', 'cache_meta');
  if ~isfield(S, 'report') || ~isfield(S, 'cache_meta')
    fprintf('[demo_isdf_scaling] stale cache ignored: %s\n', cache_path);
    return
  end
  cur = local_cache_meta(UC_save_dir, types, ratio, config, data_sc);
  if ~local_cache_meta_equal(S.cache_meta, cur)
    fprintf('[demo_isdf_scaling] stale cache ignored: %s\n', cache_path);
    return
  end
  fprintf('[demo_isdf_scaling] load cache: %s\n', cache_path);
  report = S.report;
  report.from_cache = true;
  report.cache_path = cache_path;
end

function path = local_cache_path(UC_save_dir, types, ratio)
  if numel(types) == 1
    name = sprintf('demo_isdf_scaling_%s_k%d_%d_%d.mat', ...
      types{1}, ratio(1), ratio(2), ratio(3));
  else
    name = sprintf('demo_isdf_scaling_k%d_%d_%d.mat', ratio(1), ratio(2), ratio(3));
  end
  path = fullfile(UC_save_dir, name);
end

function meta = local_cache_meta(UC_save_dir, types, ratio, config, data_sc)
  meta = struct();
  meta.version = 6;
  meta.UC_save_dir = char(string(UC_save_dir));
  meta.types = cellfun(@char, types, 'UniformOutput', false);
  meta.ratio = double(ratio(:)).';
  meta.type_ck = struct();
  for it = 1:numel(types)
    typ = types{it};
    ck = local_checkpoint(UC_save_dir, typ);
    d = dir(ck);
    meta.type_ck.(typ) = struct('file', ck, 'datenum', d(1).datenum);
  end
  meta.sc_fft = double([data_sc.sys.n1, data_sc.sys.n2, data_sc.sys.n3]);
  meta.nb = size(data_sc.psig{1}, 2);
  if isfield(config, 'SYSTEM') && isfield(config.SYSTEM, 'energy_band_index_max')
    meta.energy_band_index_max = double(config.SYSTEM.energy_band_index_max);
  else
    meta.energy_band_index_max = meta.nb;
  end
  if isfield(config, 'FORMAL') && isfield(config.FORMAL, 'wf_seed')
    meta.wf_seed = double(config.FORMAL.wf_seed);
  else
    meta.wf_seed = NaN;
  end
  meta.qwq_ratio = struct();
  for it = 1:numel(types)
    typ = types{it};
    meta.qwq_ratio.(typ) = local_target_isdf_ratio(config, typ);
  end
end

function tf = local_cache_meta_equal(a, b)
  keys = {'version', 'UC_save_dir', 'types', 'ratio', 'sc_fft', 'nb', ...
    'energy_band_index_max', 'wf_seed', 'qwq_ratio'};
  for k = 1:numel(keys)
    key = keys{k};
    if ~isfield(a, key) || ~isfield(b, key) || ~isequal(a.(key), b.(key))
      tf = false;
      return
    end
  end
  if ~isfield(a, 'type_ck') || ~isfield(b, 'type_ck')
    tf = false;
    return
  end
  fn = fieldnames(a.type_ck);
  if ~isequal(sort(fn), sort(fieldnames(b.type_ck)))
    tf = false;
    return
  end
  for k = 1:numel(fn)
    typ = fn{k};
    if ~strcmp(a.type_ck.(typ).file, b.type_ck.(typ).file) ...
        || a.type_ck.(typ).datenum ~= b.type_ck.(typ).datenum
      tf = false;
      return
    end
  end
  tf = true;
end

function ratio = local_target_isdf_ratio(config, isdf_type)
  if ~isfield(config, 'ISDF') || ~isstruct(config.ISDF)
    ratio = 0;
    return;
  end
  cfg = config.ISDF;
  switch lower(strtrim(char(string(isdf_type))))
    case 'vc'
      field = 'isdf_ratio_type1';
    case 'vn'
      field = 'isdf_ratio_type2';
    case 'nn'
      field = 'isdf_ratio_type3';
    otherwise
      field = 'isdf_ratio';
  end
  if isfield(cfg, field) && ~isempty(cfg.(field)) && double(cfg.(field)) > 0
    ratio = double(cfg.(field));
  elseif isfield(cfg, 'isdf_ratio') && double(cfg.isdf_ratio) > 0
    ratio = double(cfg.isdf_ratio);
  else
    ratio = 0;
  end
end

function fpath = local_checkpoint(source_dir, isdf_type)
  patt = fullfile(source_dir, sprintf('isdf_adaptive_checkpoint_%s_id*.mat', isdf_type));
  d = dir(patt);
  if isempty(d)
    error('demo_isdf_scaling:MissingCheckpoint', 'No %s under %s.', patt, source_dir);
  end
  [~, ord] = sort([d.datenum], 'descend');
  fpath = fullfile(d(ord(1)).folder, d(ord(1)).name);
end

function local_ensure_mex()
  need = {
    'isdf.adaptive.isdf_schur_rank1_mex', @() isdf.adaptive.build_schur_rank1_mex
    'isdf.prod_mex', @() isdf.build_prod_mex
    'isdf.adaptive.isdf_schur_rank1_prod_mex', @() isdf.adaptive.build_schur_rank1_prod_mex
  };
  for k = 1:size(need, 1)
    if isempty(which(need{k, 1}))
      fprintf('Building missing MEX: %s\n', need{k, 1});
      need{k, 2}();
    end
  end
end
