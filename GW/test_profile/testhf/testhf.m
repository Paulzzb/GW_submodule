def = filename_map();
load(def.GWinput);
GWgroundstate.psir = get_wavefunc_real(GWgroundstate.psig, GWgroundstate.Ggrid4psig);


Ex = gw_x(GWgroundstate, config);

config.ISDF.isisdf = 0;
Ex_ = gw_x(GWgroundstate, config);
[Ex, abs(Ex - Ex_) ./ abs(Ex)]
