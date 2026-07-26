% RUN_ISDF_VS_ISDFTEST  Compare +isdf and +isdftest via input_driver / isdftest.driver.
%
% Phase 1: input_driver('./test') -> service_driver -> isdf.driver
%          HF reports archived to results/isdf/
%          adaptive checkpoints remain in case SAVE/ (isdf_adaptive_checkpoint_*.mat)
% Phase 2: relay restore -> isdftest.driver (reports in results/isdftest/)
%          reuses Phase 1 adaptive checkpoints via import_isdf_adaptive_root
%
% Uses the MATLAB current directory as the case directory (must contain test).
% From any case directory in MATLAB:
%   cd /path/to/case
%   run_isdf_vs_isdftest
%
% Optional: run_isdf_vs_isdftest(true) to rebuild SAVE/ from groundstate.

function run_isdf_vs_isdftest(force_rebuild)
  cdir = pwd;
  validate_case_dir(cdir);

  gw_root = fileparts(fileparts(mfilename('fullpath')));
  cd(gw_root);
  QPstartup();
  cd(cdir);

  if nargin < 1 || isempty(force_rebuild)
    force_rebuild = false;
  end

  out_isdf = fullfile(cdir, 'results', 'isdf');
  out_isdftest = fullfile(cdir, 'results', 'isdftest');
  if force_rebuild && exist(fullfile(cdir, 'SAVE'), 'dir')
    fprintf('Removing cached SAVE/ for rebuild.\n');
    rmdir(fullfile(cdir, 'SAVE'), 's');
  end
  if exist(out_isdf, 'dir')
    rmdir(out_isdf, 's');
  end
  if exist(out_isdftest, 'dir')
    rmdir(out_isdftest, 's');
  end
  mkdir(out_isdf);
  mkdir(out_isdftest);

  service_reset_persistent();
  packages_reset_persistent();
  clear_local_hf_reports(cdir);

  fprintf('\n=== Phase 1: input_driver (+isdf via service_driver) ===\n');
  t0 = tic;
  input_driver('./test');
  fprintf('Phase 1 wall time: %.1f s\n', toc(t0));

  n_isdf = archive_hf_reports(cdir, out_isdf, '+isdf');
  fprintf('Archived %d HF/adaptive report file(s) under %s\n', n_isdf, out_isdf);

  stage_path = fullfile(cdir, 'test_relay_stage.mat');
  if ~isfile(stage_path)
    stage_path = fullfile(cdir, 'SAVE', 'test_relay_stage.mat');
  end
  if ~isfile(stage_path)
    error('run_isdf_vs_isdftest:MissingStage', ...
      'test_relay_stage.mat not found after input_driver.');
  end

  cfg_path = fullfile(cdir, 'SAVE', 'config.mat');
  if ~isfile(cfg_path)
    error('run_isdf_vs_isdftest:MissingConfig', 'SAVE/config.mat not found.');
  end
  load(cfg_path, 'config');
  if ~config.ISDF.isisdf
    error('run_isdf_vs_isdftest:ISDFOff', 'Set isisdf=1 in test input.');
  end

  fprintf('\n=== Phase 2: isdftest.driver (+isdftest) ===\n');
  ensure_isdftest_mex();
  relay.stage_from_db(stage_path);
  relay.restore();
  isdftest.debug.init_from_config(config);

  % Reuse +isdf adaptive checkpoints (written under case SAVE/ in Phase 1).
  % isdftest.adaptive.adaptiveisdf reads config.ISDF.import_isdf_adaptive_root and
  % imports isdf_adaptive_checkpoint_<desc>_id<N>.mat instead of rerunning adaptive.
  isdf_adaptive_save = fullfile(cdir, 'SAVE');
  config.ISDF.import_isdf_adaptive_root = isdf_adaptive_save;
  log_isdf_adaptive_checkpoints(isdf_adaptive_save);

  oldpwd = pwd;
  cleanup = onCleanup(@() cd(oldpwd)); %#ok<NASGU>
  cd(out_isdftest);
  t1 = tic;
  isdftest.driver([], config);
  fprintf('Phase 2 wall time: %.1f s\n', toc(t1));

  n_isdftest = numel(dir(fullfile(out_isdftest, 'isdf_validate_HF_id*.txt')));
  fprintf('Wrote %d isdf_validate_HF_id*.txt under %s\n', n_isdftest, out_isdftest);

  hist_paths = generate_diff_histograms(cdir, out_isdf, out_isdftest);
  summary_path = write_comparison_summary(cdir, out_isdf, out_isdftest);
  fprintf('\nComparison summary: %s\n', summary_path);
  if ~isempty(hist_paths)
    fprintf('Histogram outputs (%d):\n', numel(hist_paths));
    for i = 1:numel(hist_paths)
      fprintf('  %s\n', hist_paths{i});
    end
  end
  fprintf('run_isdf_vs_isdftest: done.\n');
end

function validate_case_dir(cdir)
  test_path = fullfile(cdir, 'test');
  if ~isfile(test_path)
    error('run_isdf_vs_isdftest:MissingTest', ...
      ['Current directory must contain a test input file.\n' ...
       '  pwd: %s\n' ...
       '  expected: %s'], cdir, test_path);
  end
end

function clear_local_hf_reports(cdir)
  patterns = {'isdf_validate_HF_id*.txt', 'adaptiveisdf_id*.txt'};
  for k = 1:numel(patterns)
    d = dir(fullfile(cdir, patterns{k}));
    for j = 1:numel(d)
      delete(fullfile(cdir, d(j).name));
    end
  end
end

function n = archive_hf_reports(src_dir, dst_dir, label)
  if exist(dst_dir, 'dir') ~= 7
    mkdir(dst_dir);
  end
  patterns = {'isdf_validate_HF_id*.txt', 'adaptiveisdf_id*.txt'};
  n = 0;
  for k = 1:numel(patterns)
    d = dir(fullfile(src_dir, patterns{k}));
    for j = 1:numel(d)
      src = fullfile(src_dir, d(j).name);
      dst = fullfile(dst_dir, d(j).name);
      copyfile(src, dst);
      delete(src);
      n = n + 1;
    end
  end
  stamp = fullfile(dst_dir, 'RUN_INFO.txt');
  fid = fopen(stamp, 'w');
  if fid >= 0
    fprintf(fid, 'package: %s\n', label);
    fprintf(fid, 'archived: %s\n', datestr(now, 31));
    fprintf(fid, 'source_dir: %s\n', src_dir);
    fclose(fid);
    n = n + 1;
  end
end

function fpath = write_comparison_summary(cdir, out_isdf, out_isdftest)
  fpath = fullfile(cdir, 'results', 'comparison_summary.txt');
  fid = fopen(fpath, 'w');
  if fid < 0
    error('run_isdf_vs_isdftest:Summary', 'Cannot write %s', fpath);
  end
  fprintf(fid, 'ISDF vs ISDFTEST HF validation comparison\n');
  fprintf(fid, 'Generated: %s\n\n', datestr(now, 31));

  write_report_list(fid, '+isdf (input_driver / isdf.driver)', out_isdf);
  fprintf(fid, '\n');
  write_report_list(fid, '+isdftest (isdftest.driver)', out_isdftest);
  fclose(fid);
end

function write_report_list(fid, title, dirpath)
  fprintf(fid, '--- %s ---\n', title);
  fprintf(fid, 'directory: %s\n', dirpath);
  hf = dir(fullfile(dirpath, 'isdf_validate_HF_id*.txt'));
  if isempty(hf)
    fprintf(fid, '(no isdf_validate_HF_id*.txt)\n');
    return;
  end
  for k = 1:numel(hf)
    fprintf(fid, '  %s\n', hf(k).name);
    extract_global_sums(fid, fullfile(dirpath, hf(k).name));
  end
end

function extract_global_sums(fid, fpath)
  txt = fileread(fpath);
  tok = regexp(txt, '"sum\(E_HF\)"\s+"([^"]+)"', 'tokens', 'once');
  if ~isempty(tok)
    fprintf(fid, '    sum(E_HF) = %s\n', tok{1});
  end
  tok = regexp(txt, '"sum\(E_HF_ISDF\)"\s+"([^"]+)"', 'tokens', 'once');
  if ~isempty(tok)
    fprintf(fid, '    sum(E_HF_ISDF) = %s\n', tok{1});
  end
  tok = regexp(txt, '"sum\(E_HF\) - sum\(E_HF_ISDF\)"\s+"([^"]+)"', 'tokens', 'once');
  if ~isempty(tok)
    fprintf(fid, '    sum diff = %s\n', tok{1});
  end
end

function ensure_isdftest_mex()
  local_build_if_missing('isdftest.prod_mex', @() isdftest.build_prod_mex());
  local_build_if_missing('isdftest.adaptive.isdf_schur_rank1_mex', ...
    @() isdftest.adaptive.build_schur_rank1_mex());
  local_build_if_missing('isdftest.adaptive.isdf_schur_rank1_prod_mex', ...
    @() isdftest.adaptive.build_schur_rank1_prod_mex());
end

function local_build_if_missing(symbol_name, build_fn)
  if isempty(which(symbol_name))
    fprintf('Building missing MEX: %s\n', symbol_name);
    build_fn();
  end
end

function log_isdf_adaptive_checkpoints(isdf_save_dir)
  d = dir(fullfile(isdf_save_dir, 'isdf_adaptive_checkpoint_*.mat'));
  fprintf('Phase 2 adaptive import root: %s\n', isdf_save_dir);
  if isempty(d)
    warning('run_isdf_vs_isdftest:NoIsdfAdaptiveCheckpoint', ...
      ['No isdf_adaptive_checkpoint_*.mat under %s; Phase 2 will rerun adaptiveisdf.'], ...
      isdf_save_dir);
    return;
  end
  fprintf('Phase 2 will import %d +isdf adaptive checkpoint(s):\n', numel(d));
  for k = 1:numel(d)
    fprintf('  %s\n', d(k).name);
  end
end

function out_paths = generate_diff_histograms(cdir, out_isdf, out_isdftest)
% Generate dEx histograms from HF validation reports (isdf and isdftest).
  out_paths = {};
  hist_dir = fullfile(cdir, 'results', 'hist');
  if exist(hist_dir, 'dir') ~= 7
    mkdir(hist_dir);
  end

  entries = [ ...
    collect_hist_entries('isdf', out_isdf), ...
    collect_hist_entries('isdftest', out_isdftest) ...
  ];
  if isempty(entries)
    warning('run_isdf_vs_isdftest:NoHistSource', ...
      'No isdf_validate_HF_id*.txt found; skip histogram generation.');
    return;
  end

  ts = datestr(now, 'yyyymmdd_HHMMSS');
  for i = 1:numel(entries)
    e = entries(i);
    if strcmp(e.pkg, 'isdf')
      [~, ~, ~, report] = isdf.validation.isdf_validate_HF(e.id, out_isdf);
    else
      [~, ~, ~, report] = isdftest.validation.isdf_validate_HF(e.id, out_isdftest);
    end
    dEx = report.E_diff(:);
    dEx = dEx(isfinite(dEx));
    if isempty(dEx)
      warning('run_isdf_vs_isdftest:EmptyDiff', ...
        'Empty E_diff for %s id=%d (%s); skip.', e.pkg, int32(e.id), e.desc);
      continue;
    end

    fig_h = figure('Color', 'w', 'Visible', 'off', ...
      'Name', sprintf('%s-%s-id%d histogram', e.pkg, e.desc, int32(e.id)));
    histogram(dEx, 40);
    grid on;
    title(sprintf('%s | %s (id=%d)', upper(e.pkg), e.desc, int32(e.id)));
    xlabel('dEx = E\_HF - E\_HF\_ISDF (eV)');
    ylabel('Count');

    base = sprintf('%s_%s_id%d_hist_%s', e.pkg, e.desc, int32(e.id), ts);
    png_path = fullfile(hist_dir, [base, '.png']);
    fig_path = fullfile(hist_dir, [base, '.fig']);
    try
      exportgraphics(fig_h, png_path, 'Resolution', 220);
    catch
      saveas(fig_h, png_path);
    end
    savefig(fig_h, fig_path);
    close(fig_h);
    out_paths{end + 1} = png_path; %#ok<AGROW>
    out_paths{end + 1} = fig_path; %#ok<AGROW>
  end
end

function entries = collect_hist_entries(pkg, report_dir)
% Discover id + desc from existing HF report text files.
  entries = struct('pkg', {}, 'id', {}, 'desc', {});
  d = dir(fullfile(report_dir, 'isdf_validate_HF_id*.txt'));
  for k = 1:numel(d)
    name = d(k).name;
    tok = regexp(name, 'id(\d+)\.txt$', 'tokens', 'once');
    if isempty(tok)
      continue;
    end
    id = int32(str2double(tok{1}));
    txt = fileread(fullfile(report_dir, name));
    dtok = regexp(txt, '^\s*description:\s*([A-Za-z0-9_]+)\s*$', ...
      'tokens', 'once', 'lineanchors');
    if isempty(dtok)
      desc = sprintf('id%d', int32(id));
    else
      desc = lower(strtrim(dtok{1}));
    end
    entries(end + 1) = struct('pkg', pkg, 'id', id, 'desc', desc); %#ok<AGROW>
  end
end
