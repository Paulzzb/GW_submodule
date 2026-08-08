% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/08 ZZ

%RUN_ALL  Run Si_gamma verification suites under example/cases/.
%
%   From example/ in MATLAB:
%     run_all
%
%   Suites (Si bulk / gamma; GS shared from cases/Si_gamma/qe.save):
%     G1 ISDF methods (fdep=2):     Si_gamma, Si_gamma_qrcp, Si_gamma_kmeans
%     G2 ISDF vs dense (fdep=2):    Si_gamma, Si_gamma_dir
%     G3 Cauchy on/off:
%         fdep=2:   Si_gamma_ff_cauchy, Si_gamma_ff_nocauchy
%         fdep=-2:  Si_gamma_cohsex_cauchy, Si_gamma_cohsex_nocauchy
%     Plus Si_k multi-k smoke.
%
%   Physics dispatch: name contains 'gamma' → run_cohsex; else → run_gw_x.
%   Before each case: wipe SAVE/ and prior run artifacts.
%   Collect energies afterwards with collect_qp_energies.

tests_dir = fileparts(mfilename('fullpath'));
if isempty(tests_dir)
  tests_dir = pwd;
end
gw_root = fileparts(tests_dir);
addpath(tests_dir);

% G1: ISDF methods (fdep=2) | G2: + Si_gamma_dir | G3: Cauchy ff / cohsex
cases = { ...
  'Si_gamma', ...
  'Si_gamma_qrcp', ...
  'Si_gamma_kmeans', ...
  'Si_gamma_dir', ...
  'Si_gamma_ff_cauchy', ...
  'Si_gamma_ff_nocauchy', ...
  'Si_gamma_cohsex_cauchy', ...
  'Si_gamma_cohsex_nocauchy', ...
  'Si_k' ...
  };

fprintf('=== example/run_all ===\n');
fprintf('gw_root = %s\n', gw_root);

n_ok = 0;
n_skip = 0;
n_fail = 0;

for i = 1:numel(cases)
  name = cases{i};
  case_dir = fullfile(tests_dir, 'cases', name);
  fprintf('\n----- [%d/%d] %s -----\n', i, numel(cases), name);

  if ~isfolder(case_dir)
    fprintf('FAIL: missing case dir %s\n', case_dir);
    n_fail = n_fail + 1;
    continue
  end

  clean_case_outputs(case_dir);

  is_gamma = ~isempty(regexp(name, 'gamma', 'once'));

  try
    if is_gamma
      run_cohsex(case_dir);
    else
      run_gw_x(case_dir);
    end
    fprintf('OK: %s\n', name);
    n_ok = n_ok + 1;
  catch ME
    fprintf('FAIL: %s\n%s\n', name, getReport(ME, 'basic'));
    n_fail = n_fail + 1;
  end
end

fprintf('\n=== summary: ok=%d skip=%d fail=%d ===\n', n_ok, n_skip, n_fail);
if n_fail > 0
  error('run_all:FailedCases', '%d case(s) failed.', n_fail);
end
