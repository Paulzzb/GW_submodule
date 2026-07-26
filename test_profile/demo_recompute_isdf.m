function id_new = demo_recompute_isdf(case_dir, isdf_type, output_dir, inputfile)
%DEMO_RECOMPUTE_ISDF  Re-read test -> refresh config -> recompute one ISDF type (vc/vn/nn).
%
%   demo_recompute_isdf()
%   demo_recompute_isdf(case_dir)
%   demo_recompute_isdf(case_dir, isdf_type)
%   demo_recompute_isdf(case_dir, isdf_type, output_dir)
%   demo_recompute_isdf(case_dir, isdf_type, output_dir, inputfile)
%
%   isdf_type — 'vc' | 'vn' | 'nn'  (type1 / type2 / type3)
%
% Workflow:
%   1. Parse ./test (or inputfile); defaults from SAVE/data.mat; save SAVE/config.mat.
%   2. relay.restore() for base service state (no full isdf.driver).
%   3. Recompute the requested type: gen_coeff -> adaptiveisdf -> gen_tildeVq.
%
% Prerequisites: test, SAVE/data.mat, test_relay_stage.mat (or SAVE/relay_stage.mat).
% Does NOT call isdf.free() (other ISDF slots in the same session are kept).

  if nargin < 1 || isempty(case_dir)
    case_dir = pwd;
  end
  case_dir = char(string(case_dir));
  if ~isfolder(case_dir)
    error('demo_recompute_isdf:case_dir', 'Case directory not found: %s', case_dir);
  end

  if nargin < 2 || isempty(isdf_type)
    isdf_type = 'vc';
  end
  meta = local_resolve_type(isdf_type);

  if nargin < 3 || isempty(output_dir)
    output_dir = case_dir;
  end
  output_dir = char(string(output_dir));

  if nargin < 4 || isempty(inputfile)
    inputfile = fullfile(case_dir, 'test');
  else
    inputfile = char(string(inputfile));
  end
  if ~isfile(inputfile)
    error('demo_recompute_isdf:inputfile', 'Input file not found: %s', inputfile);
  end

  profile_dir = fileparts(mfilename('fullpath'));
  gw_root = fileparts(profile_dir);
  service_dir = fullfile(gw_root, 'service');

  here_before = pwd;
  cleanup_cd = onCleanup(@() cd(here_before));

  cd(service_dir);
  addpath(genpath(service_dir));
  rehash;
  local_ensure_mex_kernels();

  cd(gw_root);
  QPstartup;

  cd(case_dir);
  fprintf('demo_recompute_isdf: case_dir = %s, type = %s (%s)\n', ...
    case_dir, meta.desc, meta.label);
  fprintf('demo_recompute_isdf: inputfile = %s\n', inputfile);

  def = filename_map();

  config = read_input_param(inputfile);
  validate_required_params(config);
  storage_dir = fullfile(case_dir, config.CONTROL.storage_dir);
  data_path = fullfile(storage_dir, def.data);
  config_path = fullfile(storage_dir, def.config);

  if ~isfile(data_path)
    error('demo_recompute_isdf:data', 'Missing %s. Run input_driver first.', data_path);
  end
  data = load(data_path, 'data').data;
  config = set_default_param_value(config, data);

  if ~isfield(config, 'ISDF') || ~config.ISDF.isisdf
    error('demo_recompute_isdf:isisdf', 'config.ISDF.isisdf must be true.');
  end
  if isfield(config.ISDF, meta.compute_field) && ~logical(config.ISDF.(meta.compute_field))
    error('demo_recompute_isdf:compute', ...
      'Set %s = .true. in test to rebuild %s.', meta.compute_field, meta.desc);
  end

  GWoptions = [];
  if isfile(config_path)
    Sc = load(config_path, 'GWoptions');
    if isfield(Sc, 'GWoptions')
      GWoptions = Sc.GWoptions;
    end
  end
  if ~exist(storage_dir, 'dir')
    mkdir(storage_dir);
  end
  save(config_path, 'GWoptions', 'config');
  fprintf('demo_recompute_isdf: updated %s from %s\n', config_path, inputfile);

  stage_path = fullfile(case_dir, 'test_relay_stage.mat');
  if ~isfile(stage_path)
    stage_path = fullfile(case_dir, def.stage);
  end
  if ~isfile(stage_path)
    stage_path = fullfile(storage_dir, def.stage);
  end
  if ~isfile(stage_path)
    error('demo_recompute_isdf:stage', ...
      'Missing relay stage. Run input_driver first.');
  end

  service_reset_persistent();
  relay.stage_from_db(stage_path);
  relay.restore();
  isdf.debug.init_from_config(config);

  cfg = config.ISDF;
  local_free_desc_slots(meta.desc);

  if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
  end

  isdf.adaptive.isdf_schur_update('clear');
  isdf.adaptive.adaptive_weight('clear');

  id_coarse = isdf.isdf_add(meta.desc);
  isdf.set_nrange(id_coarse, config.SYSTEM);
  slot_data = isdf.get(id_coarse);
  Nnrange1 = double(slot_data.Nnrange1);
  Nnrange2 = double(slot_data.Nnrange2);
  nmu_target = cfg.(meta.ratio_field) * sqrt(Nnrange1 * Nnrange2);
  slot_data.nisdf = int32(max(1, ceil(nmu_target)));
  isdf.save2mod(slot_data, id_coarse);

  isdf.coeff.gen_coeff(cfg, id_coarse);
  isdf.coeff.print_coarse_grid_report(id_coarse);

  % Force full adaptive pass: old checkpoints are keyed by pool id only and may
  % belong to another desc (e.g. vc id1 file loaded when recomputing vn on id1).
  isdf.adaptive.adaptive_checkpoint_clear(id_coarse, meta.desc);

  cleanup_out = onCleanup(@() cd(case_dir));
  cd(output_dir);
  id_new = isdf.adaptive.adaptiveisdf(id_coarse, cfg);
  cd(case_dir);
  clear cleanup_out;

  isdf.gen_tildeVq(id_new, cfg);

  slot_final = isdf.get(id_new);
  if isempty(slot_final.helperqG)
    error('demo_recompute_isdf:helperqG', ...
      '%s slot id=%d has empty helperqG after gen_tildeVq.', meta.desc, double(id_new));
  end

  rep_ad = dir(fullfile(output_dir, sprintf('adaptiveisdf_id%d.txt', id_coarse)));
  if isempty(rep_ad)
    warning('demo_recompute_isdf:noAdaptiveReport', ...
      'No adaptiveisdf_id%d.txt in output_dir=%s', double(id_coarse), output_dir);
  else
    fprintf('demo_recompute_isdf: adaptive report: %s\n', ...
      fullfile(output_dir, rep_ad(1).name));
  end

  local_print_pool_summary(meta.desc);
  fprintf(['demo_recompute_isdf: OK type=%s (coarse id=%d -> adaptive id=%d, ', ...
    'Nisdf=%d, output_dir=%s)\n'], ...
    meta.desc, double(id_coarse), double(id_new), double(slot_final.nisdf), output_dir);

end

function meta = local_resolve_type(isdf_type)
  isdf_type = lower(strtrim(char(string(isdf_type))));
  switch isdf_type
    case 'vc'
      meta = struct('desc', 'vc', 'label', 'type1', ...
        'ratio_field', 'isdf_ratio_type1', 'compute_field', 'compute_vc');
    case 'vn'
      meta = struct('desc', 'vn', 'label', 'type2', ...
        'ratio_field', 'isdf_ratio_type2', 'compute_field', 'compute_vn');
    case 'nn'
      meta = struct('desc', 'nn', 'label', 'type3', ...
        'ratio_field', 'isdf_ratio_type3', 'compute_field', 'compute_nn');
    otherwise
      error('demo_recompute_isdf:type', ...
        'isdf_type must be ''vc'', ''vn'', or ''nn'' (got ''%s'').', isdf_type);
  end
end

function local_free_desc_slots(desc)
  L = isdf.manager('list');
  for k = 1:numel(L)
    if ~L(k).empty && L(k).assigned && strcmp(char(L(k).desc), desc)
      isdf.manager('free', L(k).id);
      fprintf('demo_recompute_isdf: freed old %s slot id=%d\n', desc, double(L(k).id));
    end
  end
end

function local_print_pool_summary(desc)
  L = isdf.manager('list');
  fprintf('demo_recompute_isdf: ISDF pool after %s recompute:\n', desc);
  for k = 1:numel(L)
    if L(k).empty
      continue
    end
    fprintf('  id=%d desc=%s assigned=%d\n', ...
      double(L(k).id), char(L(k).desc), double(L(k).assigned));
  end
end

function local_ensure_mex_kernels()
  local_build_if_missing('isdf.adaptive.isdf_schur_rank1_mex', ...
    @() isdf.adaptive.build_schur_rank1_mex);
  local_build_if_missing('isdf.prod_mex', @() isdf.build_prod_mex);
  local_build_if_missing('isdf.adaptive.isdf_schur_rank1_prod_mex', ...
    @() isdf.adaptive.build_schur_rank1_prod_mex);
end

function local_build_if_missing(symbol_name, build_fn)
  if isempty(which(symbol_name))
    fprintf('demo_recompute_isdf: building MEX %s\n', symbol_name);
    build_fn();
  end
end
