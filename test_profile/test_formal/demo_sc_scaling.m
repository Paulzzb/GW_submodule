function out = demo_sc_scaling(UC_save_dir, k1k2k3, opts)
%DEMO_SC_SCALING  Load UC groundstate; synthesize formal supercell groundstate.
%
%   out = demo_sc_scaling(UC_save_dir, k1k2k3)
%   out = demo_sc_scaling(UC_save_dir, k1k2k3, opts)
%
%   opts.workers     - parpool size (default: resolve_parpool_workers())
%   opts.use_cache   - logical, default true (ISDF cache)
%   opts.use_cauchy  - logical, default false; GW COmegaC via COmegaCstar
% Step 1: Read unit-cell groundstate from UC_save_dir (data.mat, config.mat).
% Step 2: Build formal supercell groundstate via load_formal_groundstate and k1k2k3.
% Does NOT call input_driver or service_driver.

  UC_save_dir = char(string(UC_save_dir));
  if ~isfolder(UC_save_dir)
    error('demo_sc_scaling:UC_save_dir', 'UC SAVE directory not found: %s', UC_save_dir);
  end

  ratio = local_parse_k1k2k3(k1k2k3);

  test_profile_dir = fileparts(mfilename('fullpath'));
  gw_root_dir = fileparts(fileparts(test_profile_dir));

  here = pwd;
  cleanup = onCleanup(@() cd(here));

  cd(gw_root_dir);
  QPstartup;

  def = filename_map();
  data_path = fullfile(UC_save_dir, def.data);
  config_path = fullfile(UC_save_dir, def.config);
  if ~isfile(data_path)
    error('demo_sc_scaling:data', 'Missing %s', data_path);
  end
  if ~isfile(config_path)
    error('demo_sc_scaling:config', 'Missing %s', config_path);
  end

  Sdata = load(data_path, 'data');
  data_uc = Sdata.data;
  Scfg = load(config_path, 'config');
  config_uc = Scfg.config;

  case_dir = fileparts(UC_save_dir);
  groundstate_dir = local_resolve_groundstate_dir(config_uc, case_dir);

  config_sc = local_build_formal_config(config_uc, ratio, UC_save_dir);
  config_sc = local_apply_run_opts(config_sc, opts);
  data_sc = load_formal_groundstate(groundstate_dir, config_sc);

  out = struct();
  out.UC_save_dir = UC_save_dir;
  out.k1k2k3 = ratio;
  out.groundstate_dir = groundstate_dir;
  out.data_uc = data_uc;
  out.config_uc = config_uc;
  out.data_sc = data_sc;
  out.config_sc = config_sc;

  fprintf(['demo_sc_scaling: UC [%d %d %d] -> SC [%d %d %d], ', ...
    'nb %d -> %d, ng %d -> %d\n'], ...
    data_uc.sys.n1, data_uc.sys.n2, data_uc.sys.n3, ...
    data_sc.sys.n1, data_sc.sys.n2, data_sc.sys.n3, ...
    size(data_uc.psig{1}, 2), size(data_sc.psig{1}, 2), ...
    data_uc.sys.ng, data_sc.sys.ng);

  out.isdf = demo_isdf_scaling(out);
  out.gw = demo_GW_scaling(out);

end

function ratio = local_parse_k1k2k3(k1k2k3)
  ratio = double(k1k2k3(:)).';
  if numel(ratio) ~= 3
    error('demo_sc_scaling:k1k2k3', 'k1k2k3 must have three elements [k1 k2 k3].');
  end
  if any(ratio < 1) || any(mod(ratio, 1) ~= 0)
    error('demo_sc_scaling:k1k2k3', 'k1/k2/k3 must be integers >= 1.');
  end
  ratio = int32(ratio);
end

function groundstate_dir = local_resolve_groundstate_dir(config_uc, case_dir)
  if ~isfield(config_uc, 'CONTROL') || ~isfield(config_uc.CONTROL, 'groundstate_dir')
    error('demo_sc_scaling:groundstate_dir', ...
      'config.CONTROL.groundstate_dir missing in UC config.mat.');
  end
  gs = strtrim(char(string(config_uc.CONTROL.groundstate_dir)));
  if isfolder(gs)
    groundstate_dir = gs;
  else
    groundstate_dir = fullfile(case_dir, gs);
  end
  if ~isfolder(groundstate_dir)
    error('demo_sc_scaling:groundstate_dir', ...
      'Groundstate directory not found: %s', groundstate_dir);
  end
end

function config_sc = local_build_formal_config(config_uc, ratio, UC_save_dir)
  config_sc = config_uc;
  scale = double(prod(ratio));

  if ~isfield(config_sc, 'SUPERCELL') || ~isstruct(config_sc.SUPERCELL)
    config_sc.SUPERCELL = default_param_values().SUPERCELL;
  end
  config_sc.SUPERCELL.k1 = ratio(1);
  config_sc.SUPERCELL.k2 = ratio(2);
  config_sc.SUPERCELL.k3 = ratio(3);
  config_sc.SUPERCELL.is_supercell = true;
  config_sc.SUPERCELL.use_sc_isdf = true;
  config_sc.SUPERCELL.isdf_source_dir = UC_save_dir;

  if ~isfield(config_sc, 'FORMAL') || ~isstruct(config_sc.FORMAL)
    config_sc.FORMAL = default_param_values().FORMAL;
  end

  if isfield(config_sc, 'CONTROL') && isstruct(config_sc.CONTROL)
    config_sc.CONTROL.groundstate_type = 'formal';
    config_sc.CONTROL.enable_k_points = false;
  end

  if isfield(config_sc, 'SYSTEM') && isstruct(config_sc.SYSTEM)
    if isfield(config_sc.SYSTEM, 'energy_band_index_max') ...
        && config_sc.SYSTEM.energy_band_index_max > 0
      config_sc.SYSTEM.energy_band_index_max = ...
        round(double(config_sc.SYSTEM.energy_band_index_max) * scale);
    end
    if isfield(config_sc.SYSTEM, 'number_bands_in_summation') ...
        && config_sc.SYSTEM.number_bands_in_summation > 0
      config_sc.SYSTEM.number_bands_in_summation = ...
        round(double(config_sc.SYSTEM.number_bands_in_summation) * scale);
    end
  end
end

function config_sc = local_apply_run_opts(config_sc, opts)
  workers = resolve_parpool_workers();
  use_cache = true;
  use_cauchy = false;
  if nargin >= 2 && ~isempty(opts) && isstruct(opts)
    if isfield(opts, 'workers') && ~isempty(opts.workers) && double(opts.workers) > 0
      workers = round(double(opts.workers));
    end
    if isfield(opts, 'use_cache') && ~isempty(opts.use_cache)
      use_cache = logical(opts.use_cache);
    end
    if isfield(opts, 'use_cauchy') && ~isempty(opts.use_cauchy)
      use_cauchy = logical(opts.use_cauchy);
    end
  end
  if ~isfield(config_sc, 'ISDFCauchy') || ~isstruct(config_sc.ISDFCauchy)
    config_sc.ISDFCauchy = default_ISDFCauchy();
  end
  config_sc.ISDFCauchy.isCauchy = use_cauchy;
  if ~isfield(config_sc.ISDFCauchy, 'optionsCauchy') ...
      || isempty(config_sc.ISDFCauchy.optionsCauchy) ...
      || ~isstruct(config_sc.ISDFCauchy.optionsCauchy)
    config_sc.ISDFCauchy.optionsCauchy = struct( ...
      'froErr', config_sc.ISDFCauchy.froErr, ...
      'MaxIter', config_sc.ISDFCauchy.MaxIter);
  end
  if ~isfield(config_sc, 'PARALLEL') || ~isstruct(config_sc.PARALLEL)
    config_sc.PARALLEL = struct('enabled', true, 'workers', int32(workers));
  else
    config_sc.PARALLEL.enabled = true;
    config_sc.PARALLEL.workers = int32(workers);
  end
  if ~isfield(config_sc, 'TESTFUNC') || ~isstruct(config_sc.TESTFUNC)
    config_sc.TESTFUNC = struct();
  end
  config_sc.TESTFUNC.use_cache = use_cache;
end
