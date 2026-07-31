% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/07/31 ZZ

function Ex = run_gw_x(case_dir)
%RUN_GW_X  Shared launcher: input_driver + gw_x_k_packages for a case directory.
%
%   Ex = run_gw_x(CASE_DIR)
%   Ex = run_gw_x()            % CASE_DIR = pwd
%
%   CASE_DIR must contain:
%     ./test          namelist
%     groundstate as given by CONTROL.groundstate_dir (usually ./qe.save)
%
%   Writes under CASE_DIR/SAVE/ via input_driver.

if nargin < 1 || isempty(case_dir)
  case_dir = pwd;
end

if ~isfolder(case_dir)
  error('run_gw_x:MissingCaseDir', 'Case directory not found: %s', case_dir);
end

old_dir = pwd;
cd(case_dir);
case_dir = pwd;
cd(old_dir);

tests_dir = fileparts(mfilename('fullpath'));
gw_root = fileparts(tests_dir);

test_file = fullfile(case_dir, 'test');
if ~isfile(test_file)
  error('run_gw_x:MissingTest', 'Missing namelist: %s', test_file);
end

gs_rel = local_groundstate_rel_from_test(case_dir);
if ~isempty(gs_rel)
  gs_dir = fullfile(case_dir, gs_rel);
  if ~local_groundstate_ready(gs_dir)
    error('run_gw_x:MissingGroundstate', ...
      ['Groundstate not ready: %s\n', ...
       'Copy QE *.save contents into that folder (see case README).'], gs_dir);
  end
end

cd(gw_root);
QPstartup;
cd(case_dir);

service_reset_persistent();
packages_reset_persistent();

fprintf('\n[run_gw_x] case_dir = %s\n', case_dir);
fprintf('[run_gw_x] Running input_driver ...\n');
t_in = tic;
input_driver('./test');
load(fullfile(case_dir, 'SAVE', 'config.mat'), 'config');
wall_input = toc(t_in);

fprintf('[run_gw_x] Running gw_x_k_packages(config) ...\n');
t_x = tic;
Ex = gw_x_k_packages(config);
wall_x = toc(t_x);

fprintf('\n[run_gw_x] Wall clock:\n');
fprintf('  input_driver + load(config): %.3f s\n', wall_input);
fprintf('  gw_x_k_packages:             %.3f s\n', wall_x);
fprintf('  total:                       %.3f s\n', wall_input + wall_x);
fprintf('[run_gw_x] size(Ex) = [%d %d], ||Ex||_F = %.6e\n', ...
  size(Ex, 1), size(Ex, 2), norm(Ex, 'fro'));

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
  if isfile(fullfile(gs_dir, 'data-file-schema.xml'))
    ok = true;
    return
  end
  d = dir(fullfile(gs_dir, 'wfc*.hdf5'));
  ok = ~isempty(d);
end
