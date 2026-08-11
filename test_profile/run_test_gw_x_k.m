function E = run_test_gw_x_k(case_dir, profileflag, groundstate_rel)
%RUN_TEST_GW_X_K  Launcher for material GW cases: ISDF vs dense Ex via test_gw_x_k.
%
%   E = run_test_gw_x_k()
%   E = run_test_gw_x_k(case_dir)
%   E = run_test_gw_x_k(case_dir, profileflag)
%   E = run_test_gw_x_k(case_dir, profileflag, groundstate_rel)
%   E = run_test_gw_x_k(profileflag)   % case_dir = pwd (1st arg logical/scalar)
%
%   case_dir          — directory containing ./test and (after input_driver) ./SAVE
%   profileflag       — pass true to enable MATLAB profiler in test_gw_x_k
%   groundstate_rel   — optional relative path under case_dir (e.g. 'QE/GaN.save');
%                       if nonempty, must exist before the run proceeds
%
%   After test_gw_x_k (via report_ex_benchmark):
%     1. Print/write case_dir/Ex.dat — ik tables + config.ISDF (inv_strategy, order first)
%     2. Save 2x2 figure under case_dir/log/
%   E.ex_dat / E.plot_png / E.plot_fig hold output paths.

  if nargin < 3
    groundstate_rel = '';
  end
  if nargin < 2
    profileflag = false;
  end
  [case_dir, profileflag, groundstate_rel] = local_parse_args( ...
    case_dir, profileflag, groundstate_rel);

  if ~isfolder(case_dir)
    error('run_test_gw_x_k:MissingCaseDir', ...
      'Case directory not found: %s', case_dir);
  end

  test_file = fullfile(case_dir, 'test');
  if ~isfile(test_file)
    error('run_test_gw_x_k:MissingTest', ...
      'Missing test input file: %s', test_file);
  end

  if isempty(groundstate_rel)
    groundstate_rel = local_groundstate_rel_from_test(case_dir);
  end

  if ~isempty(groundstate_rel)
    gs_dir = local_resolve_groundstate_dir(case_dir, groundstate_rel);
    if ~isfolder(gs_dir)
      error('run_test_gw_x_k:groundstate', ...
        'Missing ground state: %s', gs_dir);
    end
  end

  tag = local_case_tag(case_dir);
  fprintf('[run_test_gw_x_k] case_dir=%s\n', case_dir);
  if ~isempty(groundstate_rel)
    fprintf('[run_test_gw_x_k] groundstate_rel=%s\n', groundstate_rel);
  end

  E = local_call_test_gw_x_k(case_dir, profileflag);
  E = report_ex_benchmark(E, case_dir, tag);
end

function E = local_call_test_gw_x_k(case_dir, profileflag)
  test_profile_dir = fileparts(mfilename('fullpath'));
  oldpath = path;
  cleanup = onCleanup(@() path(oldpath));
  pp = strsplit(path, pathsep);
  pp = pp(~strcmp(pp, case_dir));
  path(strjoin(pp, pathsep));
  addpath(test_profile_dir, '-begin');
  E = test_gw_x_k(case_dir, profileflag);
end

function gs_rel = local_groundstate_rel_from_test(case_dir)
  gs_rel = '';
  test_file = fullfile(case_dir, 'test');
  if ~isfile(test_file)
    return
  end
  txt = fileread(test_file);
  tok = regexp(txt, 'groundstate_dir\s*=\s*[''"]?([^,''"\s]+)[''"]?', ...
    'tokens', 'once', 'ignorecase');
  if isempty(tok)
    return
  end
  gs_rel = strtrim(tok{1});
  gs_rel = regexprep(gs_rel, '^\./', '');
end

function gs_dir = local_resolve_groundstate_dir(case_dir, groundstate_rel)
  gs_rel = char(string(groundstate_rel));
  gs_rel = regexprep(gs_rel, '^\./', '');
  gs_dir = fullfile(case_dir, gs_rel);
  if isfolder(gs_dir)
    return
  end
  % Allow Si.save vs Si.SAVE style mismatch on case-insensitive FS.
  parent = fileparts(gs_dir);
  [~, base, ext] = fileparts(gs_dir);
  if isempty(ext)
    ext = '.save';
  end
  if ~isfolder(parent)
    return
  end
  d = dir(parent);
  names = {d([d.isdir]).name};
  names = names(~ismember(names, {'.', '..'}));
  want = lower([base ext]);
  for k = 1:numel(names)
    if strcmpi(names{k}, [base ext]) || strcmpi(names{k}, base)
      gs_dir = fullfile(parent, names{k});
      return
    end
  end
end

function [case_dir, profileflag, groundstate_rel] = local_parse_args(a1, a2, a3)
  case_dir = pwd;
  profileflag = false;
  groundstate_rel = '';

  if nargin < 1 || isempty(a1)
    return
  end

  if islogical(a1) || (isnumeric(a1) && isscalar(a1) && ~(ischar(a1) || isstring(a1)))
    profileflag = logical(a1);
    return
  end

  case_dir = char(string(a1));
  if nargin >= 2 && ~isempty(a2)
    profileflag = logical(a2);
  end
  if nargin >= 3 && ~isempty(a3)
    groundstate_rel = char(string(a3));
  end
end

function tag = local_case_tag(case_dir)
  [~, tag] = fileparts(case_dir);
  if isempty(tag)
    tag = 'run_test_gw_x_k';
  end
end
