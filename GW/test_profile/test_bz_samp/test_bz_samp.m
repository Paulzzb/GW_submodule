cd ../../
QPstartup
cd test_profile/test_bz_samp

load ../TMP_FILES/GWinput.mat

GWgroundstate.psir = get_wavefunc_real(GWgroundstate.psig, GWgroundstate.Ggrid4psig);

bz_samp = GWgroundstate.bz_samp;
gvec = GWgroundstate.gvec;
symminfo = GWgroundstate.symminfo;
symmtx = symminfo.mtrx;

bmatrix = bz_samp.bmatrix;
kptbz_ca = bz_samp.kptbz * inv(bmatrix);
kpt_ca = bz_samp.kpt * inv(bmatrix);

components = gvec.components;

nibz = bz_samp.nibz;
nbz = bz_samp.nbz;
kpt = bz_samp.kpt;
kptbz = bz_samp.kptbz;

% S(kibz) = (kbz-Go) 
for ikbz = 1:nbz
  ikibz = bz_samp.kbz2kibz_ind_kbz(ikbz);
  imtx = bz_samp.kbz2kibz_ind_rotation(ikbz);
  mtx = symmtx{imtx};
  iGo = bz_samp.iGolist(ikbz);
  Go = components(iGo, :);
  kbz = kptbz_ca(ikbz,:);
  kibz = kpt_ca(ikibz,:);
  
  tmp = kibz * mtx - (kbz - Go)
  if (norm(tmp) > 1e-6)
    error();
  end
end

