%TEST_GW_X_K_LEGACY  Old Si smoke script (GWinfo + gw_x_k_packages).
%
% Prefer: run_test_gw_x_k(pwd, false) or run_test_gw_x_k(pwd, false, 'Si.save')
% which calls GW/test_profile/test_gw_x_k.m (ISDF vs dense Ex comparison).

CPATH = fileparts(mfilename('fullpath'));
cd ../../
QPstartup
cd(CPATH)

service_reset_persistent;
packages_reset_persistent;

GWinfo = GWgroundstate;
load ./SAVE/GWinput.mat
load ./SAVE/config.mat

Ex = gw_x_k_packages(GWinfo, config);
