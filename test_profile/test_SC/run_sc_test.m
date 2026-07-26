function run_sc_test(case_dir)
%RUN_SC_TEST  Run input_driver for a test_SC case (ISDF / SC_ISDF + validate_hf).
%
%   run_sc_test()              % pwd
%   run_sc_test(case_dir)

  if nargin < 1 || isempty(case_dir)
    case_dir = pwd;
  end
  case_dir = char(string(case_dir));
  if ~isfolder(case_dir)
    error('run_sc_test:case_dir', 'Case directory not found: %s', case_dir);
  end

  test_profile_dir = fileparts(mfilename('fullpath'));
  gw_root_dir = fileparts(fileparts(test_profile_dir));
  service_dir = fullfile(gw_root_dir, 'service');

  here = pwd;
  cleanup = onCleanup(@() cd(here));

  cd(service_dir);
  addpath(genpath(service_dir));
  rehash;

  ensure_mex_kernel('isdf.adaptive.isdf_schur_rank1_mex', ...
    @() isdf.adaptive.build_schur_rank1_mex);
  ensure_mex_kernel('isdf.prod_mex', @() isdf.build_prod_mex);
  ensure_mex_kernel('isdf.adaptive.isdf_schur_rank1_prod_mex', ...
    @() isdf.adaptive.build_schur_rank1_prod_mex);

  cd(gw_root_dir);
  QPstartup;

  cd(case_dir);
  if ~isfile('test')
    error('run_sc_test:test', 'Missing ./test in %s', case_dir);
  end

  service_reset_persistent();
  packages_reset_persistent();

  fprintf('run_sc_test: case_dir = %s\n', case_dir);
  t0 = tic;
  input_driver('./test');
  fprintf('run_sc_test: input_driver finished in %.3f s\n', toc(t0));

  reps = dir(fullfile(pwd, 'isdf_validate_HF_id*.txt'));
  if isempty(reps)
    warning('run_sc_test:noReport', 'No isdf_validate_HF_id*.txt in %s', case_dir);
  else
    fprintf('run_sc_test: HF validation report(s):\n');
    for k = 1:numel(reps)
      fprintf('  %s\n', fullfile(reps(k).folder, reps(k).name));
    end
  end
end

function ensure_mex_kernel(symbol_name, build_fn)
  if isempty(which(symbol_name))
    fprintf('Building missing MEX: %s\n', symbol_name);
    build_fn();
  end
end
