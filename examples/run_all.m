% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/10 ZZ

%RUN_ALL  Public demo: ISDF vs dense (Gamma Si), under examples/cases/.
%
%   From examples/ in MATLAB:
%     run_all
%
%   Cases (shared GS cases/qe.save; each case has its own ./SAVE):
%     gamma_ff_isdf      fdep=2,  isisdf=1
%     gamma_ff_dir       fdep=2,  isisdf=0
%     gamma_cohsex_isdf  fdep=-2, isisdf=1
%     gamma_cohsex_dir   fdep=-2, isisdf=0
%
%   Per-case storage so timing includes full setup (no shared SAVE reuse).

examples_dir = fileparts(mfilename('fullpath'));
if isempty(examples_dir)
  examples_dir = pwd;
end
gw_root = fileparts(examples_dir);
addpath(examples_dir);

cases = { ...
   'gamma_ff_isdf', ...
   'gamma_ff_dir', ...
   'gamma_cohsex_isdf', ...
   'gamma_cohsex_dir' ...
  };

fprintf('=== examples/run_all ===\n');
fprintf('gw_root = %s\n', gw_root);

n_ok = 0;
n_fail = 0;

for i = 1:numel(cases)
  name = cases{i};
  case_dir = fullfile(examples_dir, 'cases', name);
  fprintf('\n----- [%d/%d] %s -----\n', i, numel(cases), name);

  if ~isfolder(case_dir)
    fprintf('FAIL: missing case dir %s\n', case_dir);
    n_fail = n_fail + 1;
    continue
  end

  % Per-case SAVE / logs / qp / isdf_report
  % clean_case_outputs(case_dir);

  try
    run_cohsex(case_dir);
    fprintf('OK: %s\n', name);
    n_ok = n_ok + 1;
  catch ME
    fprintf('FAIL: %s\n%s\n', name, getReport(ME, 'basic'));
    n_fail = n_fail + 1;
  end
end

fprintf('\n=== summary: ok=%d fail=%d ===\n', n_ok, n_fail);
if n_fail > 0
  error('run_all:FailedCases', '%d case(s) failed.', n_fail);
end
