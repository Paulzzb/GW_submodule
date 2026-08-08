% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/07/27 ZZ

%RUN_ALL  Run the official Si regression cases under tests/cases/.
%
%   From repo root in MATLAB:
%     cd tests
%     run_all
%
%   Physics dispatch:
%     *gamma  → run_cohsex
%     others  → run_gw_x  (input_driver + gw.x)
%   Before each case: wipe SAVE/ and prior run artifacts so nothing is reused
%   from a previous run (keeps ./test and ./qe.save).
%   Skips a case if qe.save is not populated yet (prints SKIP).
%   SrTiO3_k is deferred (large QE output; not in the default list).

tests_dir = fileparts(mfilename('fullpath'));
if isempty(tests_dir)
  tests_dir = pwd;
end
gw_root = fileparts(tests_dir);
addpath(tests_dir);

cases = { ...
  'Si_gamma', ...
  'Si_gamma_qrcp', ...
  'Si_gamma_kmeans', ...
  'Si_k' ...
  };
% SrTiO3_k deferred: QE save too large to ship with the repo.

fprintf('=== tests/run_all ===\n');
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

  gs_dir = fullfile(case_dir, 'qe.save');

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
