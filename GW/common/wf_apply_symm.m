function wf = wf_apply_symm(GWinfo, isc)

  % Too widely used, do not push stack


  

  nspinor = GWinfo.nspinor;

  ib = isc(1);
  ik_is = isc(2);
  ik_isp = isc(4);

  is = isc(3);
  
  if (is == 1)
    wf = GWinfo.psir(:, ib, ik_is, ik_isp);
    return
  end
  
  % Then rotation is ~= I_3
  if (nspinor > 1)
    msg = fprintf("Support for sop calculation is developing\n");
    error(msg);
  end

  if (nspinor == 1)
    indrot = GWinfo.mapping.fft_r_rot{is};
    wf = GWinfo.psir(indrot, ib, ik_is, ik_isp);
  end
  
end % EOF