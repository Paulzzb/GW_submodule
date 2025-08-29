cd ../../
QPstartup
cd test_profile/test_bz_samp

load ../TMP_FILES/GWinput.mat

GWgroundstate.psir = get_wavefunc_real(GWgroundstate.psig, GWgroundstate.Ggrid4psig);

Ex = gw_x_k(GWgroundstate, config);