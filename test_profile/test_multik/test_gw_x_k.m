cd ../../
QPstartup
cd test_profile/test_multik

load ./SAVE/GWinput.mat
load ./SAVE/config.mat

GWinfo = GWgroundstate;
test_bz_samp
GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);

Ex = gw_x_k(GWinfo, config);