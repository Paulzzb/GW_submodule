function report = run_formal_bench(case_dir, profileflag)
%RUN_FORMAL_BENCH  Formal SC_ISDF + Gamma COHSEX complexity benchmark.
%
%   report = run_formal_bench(case_dir)
%   report = run_formal_bench(case_dir, true)   % MATLAB profile on
%
% Requires groundstate_type='formal', frequency_dependence=-2, SUPERCELL+SC_ISDF.
% Profiling HTML: ./profile_service/index.html, ./profile_qp_cohsex/index.html
% (profsave + bundled CSS for local browser viewing).

  if nargin < 1 || isempty(case_dir)
    case_dir = pwd;
  end
  if nargin < 2 || isempty(profileflag)
    profileflag = [];
  end
  case_dir = char(string(case_dir));
  if ~isfolder(case_dir)
    error('run_formal_bench:case_dir', 'Case directory not found: %s', case_dir);
  end

  test_profile_dir = fileparts(mfilename('fullpath'));
  gw_root_dir = fileparts(fileparts(test_profile_dir));
  service_dir = fullfile(gw_root_dir, 'service');

  here = pwd;
  cleanup = onCleanup(@() cd(here));

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
  if ~isfile('test')
    error('run_formal_bench:test', 'Missing ./test in %s', case_dir);
  end

  service_reset_persistent();
  packages_reset_persistent();

  fprintf('run_formal_bench: case_dir = %s\n', case_dir);
  report = struct();
  report.case_dir = case_dir;
  report.generated = datestr(now, 'yyyy-mm-dd HH:MM:SS');
  report.profile_enabled = local_resolve_profileflag(profileflag);

  if report.profile_enabled
    fprintf('run_formal_bench: profiling ON (service phase)\n');
    profile clear
    profile on
  end
  t0 = tic;
  input_driver('./test');
  report.wall_service_s = toc(t0);
  if report.profile_enabled
    profile off
    report.profile_service_dir = fullfile(pwd, 'profile_service');
    local_profsave_html(profile('info'), report.profile_service_dir);
  end
  fprintf('run_formal_bench: service_driver finished in %.3f s\n', report.wall_service_s);

  config_path = fullfile(pwd, 'SAVE', 'config.mat');
  if ~isfile(config_path)
    error('run_formal_bench:config', 'Missing %s after input_driver.', config_path);
  end
  S = load(config_path, 'config');
  config = S.config;

  run_qp = true;
  if isfield(config, 'FORMAL') && isfield(config.FORMAL, 'run_qp_cohsex')
    run_qp = logical(config.FORMAL.run_qp_cohsex);
  end

  if run_qp
    if config.FREQUENCY.frequency_dependence ~= -2
      error('run_formal_bench:freq', 'FORMAL benchmark expects frequency_dependence = -2.');
    end
    if report.profile_enabled
      profile clear
      profile on
    end
    t1 = tic;
    E = qp_cohsex(config);
    report.wall_qp_cohsex_s = toc(t1);
    if report.profile_enabled
      profile off
      report.profile_qp_dir = fullfile(pwd, 'profile_qp_cohsex');
      local_profsave_html(profile('info'), report.profile_qp_dir);
    end
    report.qp = E;
    fprintf('run_formal_bench: qp_cohsex finished in %.3f s\n', report.wall_qp_cohsex_s);
  else
    report.wall_qp_cohsex_s = 0;
  end

  report.wall_total_s = report.wall_service_s + report.wall_qp_cohsex_s;
  if report.profile_enabled
    report.profile_portal = local_write_profile_portal(pwd, report);
  end
  write_formal_bench_report(report, pwd);
  fprintf('run_formal_bench: wrote %s\n', fullfile(pwd, 'formal_bench_report.txt'));
end

function tf = local_resolve_profileflag(profileflag)
  if ~isempty(profileflag)
    tf = logical(profileflag);
    return
  end
  tf = false;
  if ~isfile('test')
    return
  end
  try
    cfg = read_input_param('./test');
    if isfield(cfg, 'FORMAL') && isfield(cfg.FORMAL, 'enable_profile')
      tf = logical(cfg.FORMAL.enable_profile);
    end
  catch
    tf = false;
  end
end

function write_formal_bench_report(report, out_dir)
  fpath = fullfile(out_dir, 'formal_bench_report.txt');
  fid = fopen(fpath, 'w');
  if fid < 0
    error('run_formal_bench:report', 'Cannot open %s for writing.', fpath);
  end
  c = onCleanup(@() fclose(fid));
  fprintf(fid, '=== FORMAL benchmark report ===\n');
  fprintf(fid, 'Generated: %s\n', report.generated);
  fprintf(fid, 'case_dir: %s\n', report.case_dir);
  fprintf(fid, 'profile_enabled: %d\n\n', logical(report.profile_enabled));
  fprintf(fid, 'wall_service_s   : %.6f\n', report.wall_service_s);
  fprintf(fid, 'wall_qp_cohsex_s: %.6f\n', report.wall_qp_cohsex_s);
  fprintf(fid, 'wall_total_s   : %.6f\n', report.wall_total_s);
  if isfield(report, 'profile_service_dir')
    fprintf(fid, 'profile_service  : %s\n', report.profile_service_dir);
    fprintf(fid, '  open in browser : %s\n', fullfile(report.profile_service_dir, 'index.html'));
  end
  if isfield(report, 'profile_qp_dir')
    fprintf(fid, 'profile_qp_cohsex: %s\n', report.profile_qp_dir);
    fprintf(fid, '  open in browser : %s\n', fullfile(report.profile_qp_dir, 'index.html'));
  end
  if isfield(report, 'profile_portal')
    fprintf(fid, 'profile_portal   : %s\n', report.profile_portal);
  end
  if isfield(report, 'qp') && isstruct(report.qp) && isfield(report.qp, 'Ex')
    fprintf(fid, '\nEx size: [%s]\n', num2str(size(report.qp.Ex)));
  end
  fprintf(fid, '\n=== end ===\n');
end

function local_profsave_html(s, destDir)
  if isempty(s)
    warning('run_formal_bench:profsave', ...
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
      warning('run_formal_bench:profsave', ...
        'profsave(s, dest) failed (%s); saved to default location instead.', ...
        ME.message);
    catch ME2
      throw(ME2);
    end
  end
  formal_profile_browser_fix(destDir);
  local_write_profile_summary(s, fullfile(destDir, 'profile_summary.txt'));
end

function portal = local_write_profile_portal(case_dir, report)
  portal = fullfile(case_dir, 'profile_portal.html');
  fid = fopen(portal, 'w');
  if fid < 0
    warning('run_formal_bench:portal', 'Cannot write %s', portal);
    return
  end
  c = onCleanup(@() fclose(fid));
  fprintf(fid, '<!DOCTYPE html><html><head><meta charset="UTF-8"><title>FORMAL profile portal</title>');
  fprintf(fid, '<style>body{font-family:sans-serif;margin:2em} a{font-size:1.1em}</style></head><body>');
  fprintf(fid, '<h1>FORMAL benchmark profiles</h1><ul>');
  if isfield(report, 'profile_service_dir')
    rel = local_rel_path(case_dir, report.profile_service_dir);
    fprintf(fid, '<li><a href="%s/index.html">Service phase (%.1f s)</a></li>', ...
      rel, report.wall_service_s);
  end
  if isfield(report, 'profile_qp_dir')
    rel = local_rel_path(case_dir, report.profile_qp_dir);
    fprintf(fid, '<li><a href="%s/index.html">qp_cohsex phase (%.1f s)</a></li>', ...
      rel, report.wall_qp_cohsex_s);
  end
  fprintf(fid, '</ul></body></html>');
end

function rel = local_rel_path(base_dir, target_dir)
  try
    rel = char(string(relativepath(target_dir, base_dir)));
  catch
    [~, name] = fileparts(target_dir);
    rel = name;
  end
  rel = strrep(rel, '\', '/');
end

function local_write_profile_summary(s, fpath)
  if isempty(s) || ~isfield(s, 'FunctionTable') || isempty(s.FunctionTable)
    return
  end
  T = s.FunctionTable;
  if ~ismember('TotalTime', T.Properties.VariableNames)
    return
  end
  [~, ord] = sort(T.TotalTime, 'descend');
  T = T(ord, :);
  n = min(40, height(T));
  fid = fopen(fpath, 'w');
  if fid < 0
    return
  end
  c = onCleanup(@() fclose(fid));
  fprintf(fid, '=== Profile summary (top %d by TotalTime) ===\n\n', n);
  has_self = ismember('TotalTimeMinusChildren', T.Properties.VariableNames);
  fprintf(fid, '%-8s %-8s %-8s %s\n', 'Total[s]', 'Self[s]', 'Calls', 'Function');
  fprintf(fid, '%s\n', repmat('-', 1, 72));
  for i = 1:n
    name = char(string(T.FunctionName(i)));
    if numel(name) > 48
      name = name(1:48);
    end
    self_t = 0;
    if has_self
      self_t = T.TotalTimeMinusChildren(i);
    end
    fprintf(fid, '%8.3f %8.3f %8d %s\n', T.TotalTime(i), self_t, ...
      T.NumCalls(i), name);
  end
end

function ensure_mex_kernel(symbol_name, build_fn)
  if isempty(which(symbol_name))
    fprintf('Building missing MEX: %s\n', symbol_name);
    build_fn();
  end
end
