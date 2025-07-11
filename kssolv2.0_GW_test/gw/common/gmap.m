function [ind] = gmap(gvec, syms, ngk, itran, kgq, isortc ,isorti, sys)
kg = gvec.mill(isortc(1:ngk, :),:) + kgq;
kgr = kg / syms.mtrx{itran, :};
%% Compute the address of kgr to kgrad
kadd = g2fft_index(kgr, sys);
kgrad1 = gvec.index_vec(kadd);
ind = isorti(kgrad1);
end