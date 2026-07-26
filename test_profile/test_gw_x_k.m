function E = test_gw_x_k(case_dir, profileflag)
%TEST_GW_X_K  Compare exchange Ex: ISDF (isisdf=1) vs dense (isisdf=0).
%
%   1. input_driver('./test')  — test 文件通常 isisdf=1，建好 service/ISDF
%   2. load config，gw_x_k_packages(config)  → Ex_isdf (eV)
%   3. 再次 load config，设 config.ISDF.isisdf = 0，再算 → Ex_dir (eV)
%   4. 输出 Ex_isdf - Ex_dir 及范数统计
%
%   E = test_gw_x_k(CASE_DIR)
%   E = test_gw_x_k(CASE_DIR, true)   % profiler on Ex_isdf and Ex_dir runs

  if nargin < 1 || isempty(case_dir)
    case_dir = pwd;
  end
  if nargin < 2 || isempty(profileflag)
    profileflag = false;
  end

  if ~isfolder(case_dir)
    error('test_gw_x_k:MissingCaseDir', ...
      'Case directory not found: %s', case_dir);
  end

  test_profile_dir = fileparts(mfilename('fullpath'));
  gw_root_dir = fileparts(test_profile_dir);
  service_dir = fullfile(gw_root_dir, 'service');

  if ~isfolder(service_dir)
    error('test_gw_x_k:MissingServiceDir', ...
      'Service directory not found: %s', service_dir);
  end

  cd(service_dir);
  addpath(genpath(service_dir));
  rehash;

  ensure_mex_kernel('isdf.adaptive.isdf_schur_rank1_mex', ...
    @() isdf.adaptive.build_schur_rank1_mex);
  ensure_mex_kernel('isdf.prod_mex', @() isdf.build_prod_mex);
  ensure_mex_kernel('isdf.adaptive.isdf_schur_rank1_prod_mex', ...
    @() isdf.adaptive.build_schur_rank1_prod_mex);

  cd(gw_root_dir);
  QPstartup;
  cd(case_dir);

  service_reset_persistent();
  packages_reset_persistent();

  config_path = fullfile(case_dir, 'SAVE', 'config.mat');
  tag = local_case_tag(case_dir);

  if profileflag
    profile clear
    profile on
  end
  t_in = tic;
  input_driver('./test');
  wall_input = toc(t_in);
  if profileflag
    profile off
    local_profsave_html(profile('info'), fullfile(case_dir, 'profile_input'));
  end

  S = load(config_path, 'config');
  config_isdf = S.config;
  if ~logical(config_isdf.ISDF.isisdf)
    warning('test_gw_x_k:ExpectedIsdf', ...
      '[%s] test file has ISDF.isisdf=0; ISDF branch still labeled Ex_isdf.', tag);
  end

  fprintf('\n[%s] (1) gw_x_k_packages with isisdf=%d ...\n', tag, logical(config_isdf.ISDF.isisdf));
  [Ex_isdf, wall_ex_isdf] = local_run_ex(config_isdf, case_dir, profileflag, 'profile_ex_isdf');

  fprintf('\n[%s] (2) reload config, set ISDF.isisdf=0, gw_x_k_packages ...\n', tag);
  S = load(config_path, 'config');
  config_dir = S.config;
  config_dir.ISDF.isisdf = 0;
  packages_reset_persistent();
  [Ex_dir, wall_ex_dir] = local_run_ex(config_dir, case_dir, profileflag, 'profile_ex_dir');

  dEx = Ex_isdf - Ex_dir;

  E.Ex_isdf = Ex_isdf;
  E.Ex_dir = Ex_dir;
  E.dEx = dEx;
  E.config_isdf = config_isdf;
  E.config_dir = config_dir;
  E.wall_input = wall_input;
  E.wall_ex_isdf = wall_ex_isdf;
  E.wall_ex_dir = wall_ex_dir;

  fprintf('\n[%s] Wall clock (tic/toc):\n', tag);
  fprintf('  input_driver:              %.3f s\n', wall_input);
  fprintf('  gw_x_k (isisdf=1):         %.3f s\n', wall_ex_isdf);
  fprintf('  gw_x_k (isisdf=0):         %.3f s\n', wall_ex_dir);

  fprintf('\n[%s] Ex comparison (eV), size [%d %d]:\n', ...
    tag, size(Ex_isdf, 1), size(Ex_isdf, 2));
  fprintf('  ||Ex_isdf||_F  = %.6e\n', norm(Ex_isdf(:), 'fro'));
  fprintf('  ||Ex_dir||_F   = %.6e\n', norm(Ex_dir(:), 'fro'));
  fprintf('  ||dEx||_F      = %.6e\n', norm(dEx(:), 'fro'));
  fprintf('  max|dEx|       = %.6e\n', max(abs(dEx(:))));
  if norm(Ex_isdf(:), 'fro') > 0
    fprintf('  rel||dEx||_F  = %.6e\n', norm(dEx(:), 'fro') / norm(Ex_isdf(:), 'fro'));
  end
  if ~isempty(Ex_isdf)
    fprintf('  band %d ik=1: Ex_isdf=%.6f  Ex_dir=%.6f  dEx=%.6f\n', ...
      config_isdf.SYSTEM.energy_band_index_min, ...
      Ex_isdf(1, 1), Ex_dir(1, 1), dEx(1, 1));
  end

  test_profile_dir = fileparts(mfilename('fullpath'));
  addpath(test_profile_dir, '-begin');
  E = report_ex_benchmark(E, case_dir, tag);
end

function [Ex, wall_ex] = local_run_ex(config, case_dir, profileflag, prof_subdir)
  if profileflag
    profile clear
    profile on
  end
  t_ex = tic;
  Ex_ry = gw_x_k_packages(config);
  wall_ex = toc(t_ex);
  if profileflag
    profile off
    local_profsave_html(profile('info'), fullfile(case_dir, prof_subdir));
  end

  ry2ev = constant_map().ry2ev;
  nk = size(Ex_ry, 2);
  Ex = -double(Ex_ry(:, 1:nk)) * ry2ev;
end

function tag = local_case_tag(case_dir)
  [~, tag] = fileparts(case_dir);
  if isempty(tag)
    tag = 'test_gw_x_k';
  end
end

function ensure_mex_kernel(symbol_name, build_fn)
  if isempty(which(symbol_name))
    fprintf('Building missing MEX: %s\n', symbol_name);
    build_fn();
  end
end

function local_profsave_html(s, destDir)
  if isempty(s)
    warning('test_gw_x_k:profsave', ...
      'profile(''info'') is empty; skip HTML export for %s.', destDir);
    return
  end
  if ~isfolder(destDir)
    mkdir(destDir);
  end
  try
    profsave(s, destDir);
  catch ME
    try
      profsave(s);
      warning('test_gw_x_k:profsave', ...
        'profsave(s, dest) failed (%s); saved to default location instead.', ...
        ME.message);
    catch ME2
      throw(ME2);
    end
  end
end
