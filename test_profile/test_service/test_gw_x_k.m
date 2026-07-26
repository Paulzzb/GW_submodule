cd ../../
QPstartup
cd test_profile/test_

load ./SAVE/GWinput.mat
load ./SAVE/config.mat

GWinfo = GWgroundstate;
service_reset_persistent;
packages_reset_persistent;
% GWinfo = construct_GWinfo_tmp(GWinfo);
% test_bz_samp
% GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);

Ex = gw_x_k_packages(GWinfo, config);