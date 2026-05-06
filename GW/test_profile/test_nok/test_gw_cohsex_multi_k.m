% License-Identifier: GPL
%
% Minimal smoke test: qp_cohsex(config) -> Ex + COHSEX (gw_cohsex_multi_k inside).
%
% Usage (MATLAB):
%   cd GW/test_profile/test_nok
%   test_gw_cohsex_multi_k              % no profiler (default)
%   test_gw_cohsex_multi_k(true)        % profile input + qp_cohsex separately, HTML in profile_* /
%
% When profileflag is true: two profiler runs — (1) input_driver + load config, (2) qp_cohsex.
% HTML reports: ./profile_input/ and ./profile_qp/ under this folder (see fprintf at end).
%
% Prerequisite: ./SAVE/ will be created by input_driver. If QE load fails, copy SAVE/
% and Si.save/ from test_Si (README_test_nok.txt).

function E = test_gw_cohsex_multi_k(profileflag)
  if nargin < 1 || isempty(profileflag)
    profileflag = false;
  end

  cfile = mfilename('fullpath');
  cpath = fileparts(cfile);
  file_dir = './SAVE/';

  cd ../../
  QPstartup
  cd(cpath);

  service_reset_persistent();
  packages_reset_persistent();

  % --- Phase 1: input_driver + load config (optionally under profiler) ---
  if profileflag
    profile clear
    profile on
  end
  t_in = tic;
  input_driver('./test');
  load([file_dir, 'config.mat'], 'config');
  wall_input = toc(t_in);
  if profileflag
    profile off
    local_profsave_html(profile('info'), fullfile(cpath, 'profile_input'))
  end

  % --- Phase 2: qp_cohsex (optionally under profiler) ---
  if profileflag
    profile clear
    profile on
  end
  fprintf('\n[test_gw_cohsex_multi_k] Running qp_cohsex(config) ...\n');
  t_qp = tic;
  E = qp_cohsex(config);
  wall_qp = toc(t_qp);
  if profileflag
    profile off
    local_profsave_html(profile('info'), fullfile(cpath, 'profile_qp'))
  end

  fprintf('\n[test_gw_cohsex_multi_k] Wall clock (tic/toc):\n');
  fprintf('  input_driver + load(config): %.3f s\n', wall_input);
  fprintf('  qp_cohsex (Ex + COHSEX):      %.3f s\n', wall_qp);
  fprintf('  total:                        %.3f s\n', wall_input + wall_qp);
  if profileflag
    fprintf('[test_gw_cohsex_multi_k] Profiler HTML (open index.html in each folder):\n');
    fprintf('  %s\n', fullfile(cpath, 'profile_input'));
    fprintf('  %s\n', fullfile(cpath, 'profile_qp'));
  end

  fprintf('\n[test_gw_cohsex_multi_k] size(E.Eqp)    = [%d %d]\n', size(E.Eqp, 1), size(E.Eqp, 2));
  fprintf('[test_gw_cohsex_multi_k] size(E.Ex)     = [%d %d]\n', size(E.Ex, 1), size(E.Ex, 2));
  fprintf('[test_gw_cohsex_multi_k] size(E.Esx_x)  = [%d %d]\n', size(E.Esx_x, 1), size(E.Esx_x, 2));
  fprintf('[test_gw_cohsex_multi_k] size(E.Ecoh)   = [%d %d]\n', size(E.Ecoh, 1), size(E.Ecoh, 2));
  fprintf('[test_gw_cohsex_multi_k] ||E.Eqp||_F    = %.6e\n', norm(E.Eqp, 'fro'));
  fprintf('[test_gw_cohsex_multi_k] ||E.Esx_x||_F  = %.6e\n', norm(E.Esx_x, 'fro'));
  fprintf('[test_gw_cohsex_multi_k] ||E.Ecoh||_F   = %.6e\n', norm(E.Ecoh, 'fro'));
end

function local_profsave_html(s, destDir)
% Write profiler HTML; destDir is a folder name (this MATLAB expects struct from profile('info')).

  if isempty(s)
    warning('test_gw_cohsex_multi_k:profsave', ...
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
      warning('test_gw_cohsex_multi_k:profsave', ...
        'profsave(s, dest) failed (%s); saved to default location instead.', ME.message);
    catch ME2
      throw(ME2);
    end
  end
end
