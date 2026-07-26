% License-Identifier: GPL
%
% Minimal smoke test for GW/packages/gw_cohsex_multi_k.m
%
% Usage (under MATLAB):
%   cd GW/test_profile/test_Si
%   test_gw_cohsex_multi_k

function E = test_gw_cohsex_multi_k()
  cfile = mfilename('fullpath');
  cpath = fileparts(cfile);
  file_dir = './SAVE/';

  cd ../../
  QPstartup
  cd(cpath);

  service_reset_persistent();
  packages_reset_persistent();

  % Rebuild/refresh input and config in current session.
  input_driver('./test');
  % load([file_dir, 'GWinput.mat'], 'GWgroundstate');
  load([file_dir, 'config.mat'], 'config');


  % Keep these variables in workspace in case downstream routines
  % inspect caller/base context.
  % GWinfo = GWgroundstate; %#ok<NASGU>
  config = config; %#ok<NASGU>

  fprintf('\n[test_gw_cohsex_multi_k] Running gw_cohsex_multi_k() ...\n');
  t0 = tic;
  % [Esx_x, Ecoh] = gw_cohsex_multi_k(config);
  E = qp_cohsex(config);
  elapsed = toc(t0);

  fprintf('[test_gw_cohsex_multi_k] Done in %.3f s\n', elapsed);
  fprintf('[test_gw_cohsex_multi_k] size(Esx_x) = [%d %d]\n', size(Esx_x, 1), size(Esx_x, 2));
  fprintf('[test_gw_cohsex_multi_k] size(Ecoh)  = [%d %d]\n', size(Ecoh, 1), size(Ecoh, 2));
  fprintf('[test_gw_cohsex_multi_k] ||Esx_x||_F = %.6e\n', norm(Esx_x, 'fro'));
  fprintf('[test_gw_cohsex_multi_k] ||Ecoh||_F  = %.6e\n', norm(Ecoh, 'fro'));
end
