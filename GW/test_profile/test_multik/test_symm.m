cd ../../
QPstartup
cd test_profile/test_multik

load ./SAVE/GWinput.mat
load ./SAVE/config.mat

% GWgroundstate.psir = get_wavefunc_real(GWgroundstate.psig, GWgroundstate.Ggrid4psig);

bz_samp = GWgroundstate.bz_samp;
gvec = GWgroundstate.gvec;
symminfo = GWgroundstate.symminfo;
symmtx = symminfo.mtrx;

bmatrix = bz_samp.bmatrix;
kptbz_rlu = bz_samp.kptbz * inv(bmatrix);
kpt_rlu = bz_samp.kpt * inv(bmatrix);

components = gvec.components;

nibz = bz_samp.nibz;
nbz = bz_samp.nbz;
kpt_ca = bz_samp.kpt;
kptbz_ca = bz_samp.kptbz;

gvec = GWgroundstate.gvec;
mapping = GWgroundstate.mapping;

fftgrids = gvec.fftgrids;

%% test Rgrid rotation

for isym = 
