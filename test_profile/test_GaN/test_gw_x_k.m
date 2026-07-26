CPATH = fileparts(mfilename('fullpath'));
cd ../../
QPstartup
cd(CPATH)

service_reset_persistent;
packages_reset_persistent;
%
GWinfo = GWgroundstate;
% GWinfo = construct_GWinfo_tmp(GWinfo);
% test_bz_samp
% GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);
load ./SAVE/GWinput.mat
load ./SAVE/config.mat

Ex = gw_x_k_packages(GWinfo, config);