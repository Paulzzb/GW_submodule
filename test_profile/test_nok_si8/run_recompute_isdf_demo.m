%RUN_RECOMPUTE_ISDF_DEMO  Smoke test demo_recompute_isdf for vc, vn, nn on test_nok_si8.
%
%   cd GW/test_profile/test_nok_si8
%   run_recompute_isdf_demo
%
% Edit ./test before running to change ISDF parameters. Requires SAVE/data.mat
% and test_relay_stage.mat from a prior input_driver run.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

types = {'vc', 'vn', 'nn'};
for k = 1:numel(types)
  fprintf('\n=== run_recompute_isdf_demo: type=%s ===\n', types{k});
  t0 = tic;
  demo_recompute_isdf(here, types{k});
  fprintf('=== type=%s finished in %.3f s ===\n\n', types{k}, toc(t0));
end
fprintf('run_recompute_isdf_demo: all types OK\n');
