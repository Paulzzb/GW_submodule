function iout = findvector(kk, gvec)
%  iout = mill2nl_findvector(kk, gvec.fftgrid(1), gvec.fftgrid(2), gvec.fftgrid(3));
  iout = 1 + (kk(1) + (kk(1) < 0) * gvec.fftgrid(1) ) ...
        + (kk(2) + (kk(2) < 0) * gvec.fftgrid(2)) * gvec.fftgrid(1) ... 
        + (kk(3) + (kk(3) < 0) * gvec.fftgrid(3)) * gvec.fftgrid(1) * gvec.fftgrid(2); 
  if (iout >= 1 && iout <= gvec.nfftgridpts)
    iout = gvec.index_vec(iout);
    if iout >= 1 && iout <= gvec.ng
      if (any(kk ~= gvec.components(iout, :)))
        iout = 0;
      end
    else
        iout = 0;
    end
  else
    iout = 0;
  end
  
  return
end

