% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function run_sc(which)
%RUN_SC  Small SC_ISDF demo: Si2_uc → Si8_sc (from test_profile/test_SC).
%
%   From example/ in MATLAB:
%     run_sc           % unit cell then supercell + qp.launcher
%     run_sc('uc')     % Si2_uc only (input_driver + HF validate)
%     run_sc('sc')     % Si8_sc only (needs Si2_uc/SAVE)
%
%   Layout:
%     cases/Si2_uc/   real QE primitive cell, adaptive ISDF → SAVE
%     cases/Si8_sc/   formal 2×2×1 + SC_ISDF from Si2_uc/SAVE

  if nargin < 1 || isempty(which)
    which = 'all';
  end
  which = lower(strtrim(char(string(which))));

  example_dir = fileparts(mfilename('fullpath'));
  if isempty(example_dir)
    example_dir = pwd;
  end
  gw_root = fileparts(example_dir);
  addpath(example_dir);

  uc_dir = fullfile(example_dir, 'cases', 'Si2_uc');
  sc_dir = fullfile(example_dir, 'cases', 'Si8_sc');

  switch which
    case {'all', 'both'}
      local_run_uc(gw_root, uc_dir);
      local_run_sc(gw_root, sc_dir, uc_dir);
    case {'uc', 'si2', 'si2_uc'}
      local_run_uc(gw_root, uc_dir);
    case {'sc', 'si8', 'si8_sc'}
      local_run_sc(gw_root, sc_dir, uc_dir);
    otherwise
      error('run_sc:which', 'Unknown target ''%s'' (use all|uc|sc).', which);
  end
end

function local_run_uc(gw_root, uc_dir)
  fprintf('\n===== [run_sc] Si2_uc (unit-cell ISDF) =====\n');
  local_require_case(uc_dir);
  clean_case_outputs(uc_dir);

  here = pwd;
  cleanup = onCleanup(@() cd(here));
  cd(gw_root);
  QPstartup;
  cd(uc_dir);
  service_reset_persistent();
  packages_reset_persistent();

  t0 = tic;
  input_driver('./test');
  fprintf('[run_sc] Si2_uc input_driver done in %.2f s\n', toc(t0));
  local_list_hf_reports(uc_dir);
end

function local_run_sc(gw_root, sc_dir, uc_dir)
  fprintf('\n===== [run_sc] Si8_sc (SC_ISDF + formal COHSEX) =====\n');
  local_require_case(sc_dir);

  src = fullfile(uc_dir, 'SAVE');
  if ~isfolder(src)
    error('run_sc:MissingUcSave', ...
      'Missing %s — run run_sc(''uc'') first.', src);
  end

  % Do not wipe Si2_uc/SAVE; only clear SC outputs.
  clean_case_outputs(sc_dir);

  here = pwd;
  cleanup = onCleanup(@() cd(here));
  cd(gw_root);
  QPstartup;
  cd(sc_dir);
  service_reset_persistent();
  packages_reset_persistent();

  t0 = tic;
  input_driver('./test');
  wall_in = toc(t0);

  load(fullfile(sc_dir, 'SAVE', 'config.mat'), 'config');
  t1 = tic;
  E = qp.launcher(config);
  wall_qp = toc(t1);

  fprintf('[run_sc] Si8_sc input_driver %.2f s, qp.launcher %.2f s\n', ...
    wall_in, wall_qp);
  if isstruct(E) && isfield(E, 'Eqp')
    fprintf('[run_sc] size(E.Eqp)=[%d %d], ||Eqp||_F=%.6e\n', ...
      size(E.Eqp, 1), size(E.Eqp, 2), norm(E.Eqp, 'fro'));
  end
end

function local_require_case(case_dir)
  if ~isfolder(case_dir)
    error('run_sc:MissingCase', 'Case directory not found: %s', case_dir);
  end
  if ~isfile(fullfile(case_dir, 'test'))
    error('run_sc:MissingTest', 'Missing namelist: %s/test', case_dir);
  end
end

function local_list_hf_reports(case_dir)
  rep_dir = fullfile(case_dir, 'isdf_report');
  if ~isfolder(rep_dir)
    hits = dir(fullfile(case_dir, 'o-ISDF_HF_id*'));
  else
    hits = dir(fullfile(rep_dir, 'o-ISDF_HF_id*'));
  end
  if isempty(hits)
    fprintf('[run_sc] (no HF report files found yet)\n');
    return
  end
  fprintf('[run_sc] HF reports:\n');
  for k = 1:numel(hits)
    fprintf('  %s\n', fullfile(hits(k).folder, hits(k).name));
  end
end
