% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/10 ZZ

%RUN_ALL  Run Si_gamma verification suites under example/cases/.
%
%   From example/ in MATLAB:
%     run_all
%
%   Order matters for SAVE / ISDF checkpoint reuse:
%     1) Si_gamma hub (builds SAVE + adaptive checkpoints 16/16/16)
%     2) shared-SAVE cauchy cluster (full SAVE symlink)
%     3) Si_gamma_162416_exact first (builds vn @ ratio 24; vc/nn linked from hub)
%     4) other 162416* (vc/nn from hub, vn from 162416_exact)
%     5) 161616* (all adaptive checkpoints from hub)
%     6) independent: qrcp / kmeans / dir / Si_k
%
%   See link_shared_save.sh, link_isdf_checkpoints.sh, example_pipe.md.

tests_dir = fileparts(mfilename('fullpath'));
if isempty(tests_dir)
  tests_dir = pwd;
end
gw_root = fileparts(tests_dir);
addpath(tests_dir);

cases = { ...
  'Si_gamma', ...
  'Si_gamma_ff_cauchy', ...
  'Si_gamma_ff_nocauchy', ...
  'Si_gamma_cohsex_cauchy', ...
  'Si_gamma_cohsex_nocauchy', ...
  'Si_gamma_162416_exact', ...
  'Si_gamma_162416_sum', ...
  'Si_gamma_162416_ff', ...
  'Si_gamma_161616_exact', ...
  'Si_gamma_161616_sum', ...
  'Si_gamma_161616_ff', ...
  'Si_gamma_qrcp', ...
  'Si_gamma_kmeans', ...
  'Si_gamma_dir', ...
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
  local_prepare_links(tests_dir, name);

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

function local_prepare_links(tests_dir, name)
% Re-apply checkpoint / SAVE links after clean (clean wipes real SAVE dirs).
  link_ckpt = fullfile(tests_dir, 'link_isdf_checkpoints.sh');
  if exist(link_ckpt, 'file') ~= 2
    return
  end
  ratio_cases = { ...
    'Si_gamma_162416_exact', 'Si_gamma_162416_sum', 'Si_gamma_162416_ff', ...
    'Si_gamma_161616_exact', 'Si_gamma_161616_sum', 'Si_gamma_161616_ff' ...
    };
  if ~any(strcmp(name, ratio_cases))
    return
  end
  cmd = sprintf('bash %s --case %s', local_shell_quote(link_ckpt), local_shell_quote(name));
  [st, out] = system(cmd);
  if st ~= 0
    warning('run_all:LinkCheckpoints', 'link_isdf_checkpoints failed for %s:\n%s', name, out);
  else
    fprintf('%s', out);
  end
end

function s = local_shell_quote(p)
  s = ['''' strrep(p, '''', '''\'''''') ''''];
end
