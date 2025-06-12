def = filename_map();
load(def.GWinput);

GWgroundstate.psir = get_wavefunc_real(GWgroundstate.psig, GWgroundstate.Ggrid4psig);

type = 'ss';
nb = size(GWgroundstate.psig, 2);
nlist = 1:nb; mlist = 1:nb;
gvec = GWgroundstate.gvec; vol = GWgroundstate.vol;
optionsISDF = GWoptions.ISDFCauchy;
[ind_mu, zeta_mu] = isdf_main(type, GWgroundstate.psir, nlist, mlist, gvec, vol, optionsISDF);

