function rhog = mtxel(GWinfo, param)
% mtxel_sigma -> Calculate <nn | exp(iGr) | mm>, where mm in sum_range
%     use GWinfo.psir, GWinfo.gvec, and GWinfo.vol
% if sum_range is not specified, sum_range = size(GWinfo.psir, 2)
iGo = param.qs(1);
irot = param.qs(3);
vol = GWinfo.vol;
gvec = GWinfo.gvec;
nspinor = GWinfo.nspinor;


wf_left = wf_apply_symm(GWinfo, param.is);
wf_right = wf_apply_symm(GWinfo, param.os);



fftbox = conj(wf_left) .* wf_right;



if (nspinor == 2)
  msg = fprintf("Support for sop calculation is developing\n");
  error(msg);
end

fftbox = reshape(fftbox, GWinfo.gvec.fftgrid);
fftbox = vol * do_FFT(fftbox, gvec.fftgrid, 1);


if (irot == 1) % no rotation
  rhog = fftbox(GWinfo.mapping.Gomapping{iGo});
else
  rhog = fftbox(GWinfo.mapping.Gomapping{iGo}(GWinfo.mapping.g_rot{irot}));
end


end % EOF 
