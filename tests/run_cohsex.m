% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/07/27 ZZ

function E = run_cohsex(case_dir)
%RUN_COHSEX  Shared launcher: input_driver + qp_cohsex for a case directory.
%
%   E = run_cohsex(CASE_DIR)
%   E = run_cohsex()            % CASE_DIR = pwd
%
%   CASE_DIR must contain:
%     ./test          namelist
%     groundstate as given by CONTROL.groundstate_dir (usually ./qe.save)
%
%   Writes under CASE_DIR/SAVE/ via input_driver; qp outputs follow package defaults.

if nargin < 1 || isempty(case_dir)
  case_dir = pwd;
end

if ~isfolder(case_dir)
  error('run_cohsex:MissingCaseDir', 'Case directory not found: %s', case_dir);
end

old_dir = pwd;
cd(case_dir);
case_dir = pwd;
cd(old_dir);

tests_dir = fileparts(mfilename('fullpath'));
gw_root = fileparts(tests_dir);

test_file = fullfile(case_dir, 'test');
if ~isfile(test_file)
  error('run_cohsex:MissingTest', 'Missing namelist: %s', test_file);
end

gs_rel = local_groundstate_rel_from_test(case_dir);
if ~isempty(gs_rel)
  gs_dir = fullfile(case_dir, gs_rel);
  if ~local_groundstate_ready(gs_dir)
    error('run_cohsex:MissingGroundstate', ...
      ['Groundstate not ready: %s\n', ...
       'Copy QE *.save contents into that folder (see case README).'], gs_dir);
  end
end

cd(gw_root);
QPstartup;
cd(case_dir);

service_reset_persistent();
packages_reset_persistent();

fprintf('\n[run_cohsex] case_dir = %s\n', case_dir);
fprintf('[run_cohsex] Running input_driver ...\n');
t_in = tic;
input_driver('./test');
load(fullfile(case_dir, 'SAVE', 'config.mat'), 'config');
wall_input = toc(t_in);

fprintf('[run_cohsex] Running qp_cohsex(config) ...\n');
t_qp = tic;
E = qp_cohsex(config);
wall_qp = toc(t_qp);

fprintf('\n[run_cohsex] Wall clock:\n');
fprintf('  input_driver + load(config): %.3f s\n', wall_input);
fprintf('  qp_cohsex:                   %.3f s\n', wall_qp);
fprintf('  total:                       %.3f s\n', wall_input + wall_qp);

if isfield(E, 'Eqp')
  fprintf('[run_cohsex] size(E.Eqp) = [%d %d], ||Eqp||_F = %.6e\n', ...
    size(E.Eqp, 1), size(E.Eqp, 2), norm(E.Eqp, 'fro'));
end

end

function gs_rel = local_groundstate_rel_from_test(case_dir)
  gs_rel = '';
  txt = fileread(fullfile(case_dir, 'test'));
  tok = regexp(txt, 'groundstate_dir\s*=\s*[''"]?([^,''"\s]+)[''"]?', ...
    'tokens', 'once', 'ignorecase');
  if isempty(tok)
    return
  end
  gs_rel = regexprep(strtrim(tok{1}), '^\./', '');
end

function ok = local_groundstate_ready(gs_dir)
  ok = false;
  if ~isfolder(gs_dir)
    return
  end
  % Minimal QE-save fingerprint: data-file-schema.xml or any wfc*.hdf5
  if isfile(fullfile(gs_dir, 'data-file-schema.xml'))
    ok = true;
    return
  end
  d = dir(fullfile(gs_dir, 'wfc*.hdf5'));
  ok = ~isempty(d);
end
