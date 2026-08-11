function report = run_demo_sc_scaling_profile(UC_save_dir, k1k2k3, out_dir, profile_opts)
%RUN_DEMO_SC_SCALING_PROFILE  demo_sc_scaling with MATLAB profile + parpool.
%
%   report = run_demo_sc_scaling_profile(UC_save_dir, [4 4 4], out_dir)
%   report = run_demo_sc_scaling_profile(..., out_dir, profile_opts)
%
%   profile_opts.use_cauchy     - GW COmegaC via COmegaCstar (default false)
%   profile_opts.use_isdf_cache - reuse ISDF cache; skip delete (default false)
%
%   Workers: resolve_parpool_workers() (Slurm / host CPUs); override with PARPOOL_WORKERS.
%   Writes profile HTML to out_dir/profile_demo_sc/ and demo_sc_scaling_report.txt.

  if nargin < 1 || isempty(UC_save_dir)
    error('run_demo_sc_scaling_profile:UC', 'UC_save_dir is required.');
  end
  if nargin < 2 || isempty(k1k2k3)
    k1k2k3 = [4 4 4];
  end
  if nargin < 3 || isempty(out_dir)
    out_dir = fullfile(fileparts(mfilename('fullpath')), ...
      sprintf('LiH_%d%d%d_demo_scaling', k1k2k3(1), k1k2k3(2), k1k2k3(3)));
  end
  out_dir = char(string(out_dir));
  if ~isfolder(out_dir)
    mkdir(out_dir);
  end

  workers = local_resolve_workers();
  ratio = k1k2k3(:).';
  UC_save_dir = char(string(UC_save_dir));

  use_cauchy = false;
  use_isdf_cache = false;
  if nargin >= 4 && ~isempty(profile_opts) && isstruct(profile_opts)
    if isfield(profile_opts, 'use_cauchy') && ~isempty(profile_opts.use_cauchy)
      use_cauchy = logical(profile_opts.use_cauchy);
    end
    if isfield(profile_opts, 'use_isdf_cache') && ~isempty(profile_opts.use_isdf_cache)
      use_isdf_cache = logical(profile_opts.use_isdf_cache);
    end
  end

  cache_path = fullfile(UC_save_dir, sprintf('demo_isdf_scaling_k%d_%d_%d.mat', ...
    ratio(1), ratio(2), ratio(3)));
  if ~use_isdf_cache && isfile(cache_path)
    fprintf('[run_demo_sc_scaling_profile] remove stale cache: %s\n', cache_path);
    delete(cache_path);
  elseif use_isdf_cache && isfile(cache_path)
    fprintf('[run_demo_sc_scaling_profile] reuse ISDF cache: %s\n', cache_path);
  end

  opts = struct('workers', workers, 'use_cache', use_isdf_cache, ...
    'use_cauchy', use_cauchy);

  fprintf(['[run_demo_sc_scaling_profile] UC=%s k=[%d %d %d] workers=%d ', ...
    'out=%s use_cauchy=%d use_isdf_cache=%d\n'], ...
    UC_save_dir, ratio(1), ratio(2), ratio(3), workers, out_dir, ...
    use_cauchy, use_isdf_cache);

  local_ensure_parpool(workers);

  profile clear
  profile on
  t0 = tic;
  out = demo_sc_scaling(UC_save_dir, ratio, opts);
  wall_s = toc(t0);
  prof = profile('info');
  profile off

  profile_dir = fullfile(out_dir, 'profile_demo_sc');
  local_profsave_html(prof, profile_dir);

  report = struct();
  report.generated = datestr(now, 'yyyy-mm-dd HH:MM:SS');
  report.UC_save_dir = UC_save_dir;
  report.k1k2k3 = double(ratio);
  report.workers = workers;
  report.use_cauchy = use_cauchy;
  report.use_isdf_cache = use_isdf_cache;
  report.wall_s = wall_s;
  report.profile_dir = profile_dir;
  report.out_dir = out_dir;
  if isfield(out, 'gw') && isfield(out.gw, 'wall_s')
    report.gw_wall_s = out.gw.wall_s;
  end
  if isfield(out, 'isdf')
    report.isdf_from_cache = isfield(out.isdf, 'from_cache') && out.isdf.from_cache;
    if isfield(out.isdf, 'vc')
      report.Nmu_vc = size(out.isdf.vc.helperqG, 2);
    end
    if isfield(out.isdf, 'nn')
      report.Nmu_nn = size(out.isdf.nn.helperqG, 2);
    end
  end

  rep_path = fullfile(out_dir, 'demo_sc_scaling_report.txt');
  local_write_report(rep_path, report, prof);
  fprintf('[run_demo_sc_scaling_profile] wall=%.3f s report=%s\n', wall_s, rep_path);
  fprintf('[run_demo_sc_scaling_profile] profile=%s/index.html\n', profile_dir);
end

function n = local_resolve_workers()
  n = resolve_parpool_workers();
  fprintf('[run_demo_sc_scaling_profile] resolve_parpool_workers -> %d\n', n);
end

function local_profsave_html(s, destDir)
  if isempty(s)
    warning('run_demo_sc_scaling_profile:profsave', 'profile info empty.');
    return;
  end
  if ~isfolder(destDir)
    mkdir(destDir);
  end
  profsave(s, destDir);
  formal_profile_browser_fix(destDir);
  local_write_profile_summary(s, fullfile(destDir, 'profile_summary.txt'));
end

function local_ensure_parpool(n)
  n = max(1, round(double(n)));
  c = parcluster('local');
  if c.NumWorkers < n
    c.NumWorkers = n;
  end
  pool = gcp('nocreate');
  if isempty(pool) || pool.NumWorkers ~= n
    if ~isempty(pool)
      delete(pool);
    end
    parpool(c, n);
  end
  fprintf('[run_demo_sc_scaling_profile] parpool NumWorkers=%d\n', gcp().NumWorkers);
end

function local_write_profile_summary(s, fpath)
  if isempty(s) || ~isfield(s, 'FunctionTable') || isempty(s.FunctionTable)
    return;
  end
  T = s.FunctionTable;
  if istable(T)
    if ~ismember('TotalTime', T.Properties.VariableNames)
      return;
    end
    [~, ord] = sort(T.TotalTime, 'descend');
    T = T(ord, :);
    n = min(40, height(T));
    fid = fopen(fpath, 'w');
    if fid < 0
      return;
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
    return;
  end
  if isstruct(T) && isfield(T, 'TotalTime')
    [~, ord] = sort([T.TotalTime], 'descend');
    T = T(ord);
    n = min(40, numel(T));
    fid = fopen(fpath, 'w');
    if fid < 0
      return;
    end
    c = onCleanup(@() fclose(fid));
    fprintf(fid, '=== Profile summary (top %d by TotalTime) ===\n\n', n);
    fprintf(fid, '%-8s %-8s %-8s %s\n', 'Total[s]', 'Self[s]', 'Calls', 'Function');
    fprintf(fid, '%s\n', repmat('-', 1, 72));
    for i = 1:n
      name = char(string(T(i).FunctionName));
      if numel(name) > 48
        name = name(1:48);
      end
      self_t = 0;
      if isfield(T, 'TotalTimeMinusChildren')
        self_t = T(i).TotalTimeMinusChildren;
      end
      fprintf(fid, '%8.3f %8.3f %8d %s\n', T(i).TotalTime, self_t, ...
        T(i).NumCalls, name);
    end
  end
end

function local_write_report(fpath, report, prof)
  fid = fopen(fpath, 'w');
  if fid < 0
    warning('run_demo_sc_scaling_profile:report', 'Cannot write %s', fpath);
    return;
  end
  c = onCleanup(@() fclose(fid));
  fprintf(fid, '=== demo_sc_scaling profile run ===\n');
  fprintf(fid, 'generated   : %s\n', report.generated);
  fprintf(fid, 'UC_save_dir : %s\n', report.UC_save_dir);
  fprintf(fid, 'k1k2k3      : [%g %g %g]\n', report.k1k2k3);
  fprintf(fid, 'workers     : %d\n', report.workers);
  if isfield(report, 'use_cauchy')
    fprintf(fid, 'use_cauchy  : %d\n', report.use_cauchy);
  end
  if isfield(report, 'use_isdf_cache')
    fprintf(fid, 'use_isdf_cache : %d\n', report.use_isdf_cache);
  end
  fprintf(fid, 'wall_s      : %.6f\n', report.wall_s);
  if isfield(report, 'gw_wall_s')
    fprintf(fid, 'gw_wall_s   : %.6f\n', report.gw_wall_s);
  end
  if isfield(report, 'Nmu_vc')
    fprintf(fid, 'Nmu_vc      : %d\n', report.Nmu_vc);
  end
  if isfield(report, 'Nmu_nn')
    fprintf(fid, 'Nmu_nn      : %d\n', report.Nmu_nn);
  end
  fprintf(fid, 'profile_dir : %s\n', report.profile_dir);
  fprintf(fid, 'profile_url : %s/index.html\n', report.profile_dir);
  if ~isempty(prof) && isfield(prof, 'FunctionTable') && ~isempty(prof.FunctionTable)
    T = prof.FunctionTable;
    if istable(T) && ismember('TotalTime', T.Properties.VariableNames)
      [mx, ix] = max(T.TotalTime);
      fprintf(fid, 'profile_top : %s (%.3f s)\n', char(string(T.FunctionName(ix))), mx);
    elseif isstruct(T) && isfield(T, 'TotalTime')
      [mx, ix] = max([T.TotalTime]);
      fprintf(fid, 'profile_top : %s (%.3f s)\n', char(string(T(ix).FunctionName)), mx);
    end
  end
  fprintf(fid, '\n=== end ===\n');
end
