%RUN_RECOMPUTE_ISDF_VC_DEMO  Run demo_recompute_isdf_vc on test_nok_si8.
%
%   cd GW/test_profile/test_nok_si8
%   run_recompute_isdf_vc_demo
%
% Requires SAVE/config.mat and test_relay_stage.mat from a prior input_driver run.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
demo_recompute_isdf_vc(here);
