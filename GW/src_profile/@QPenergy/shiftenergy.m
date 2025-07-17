function obj = shiftenergy(obj)
  n = length(obj.ev);
  TOL_DEGENERACY = obj.TOL_DEGENERACY;
  freq_dep = obj.freq_dep;

  if ((freq_dep == 0) || (freq_dep == 1)) % Use Esex_x and Ecoh
    Etmp1 = obj.Esex_x;
    Etmp2 = obj.Ecoh;
  elseif (freq_dep == 2) % Use Eres and Eint
    Etmp1 = obj.Eres;
    Etmp2 = obj.Eint;
  end
  Ex = obj.Ex;


  if (size(obj.ev) ~= n)
    msg = fprintf("For multi-kpoints, still under construction!\n");
    GWerror(msg);
  end
  ndeg = zeros(n, 1);
  nl = 1;
  ndeg(nl) = 1;
  for i = 2:n
    dek = obj.ev(i, :) - obj.ev(i-1, :);
    iflag = 0; if (abs(dek) < TOL_DEGENERACY); iflag=1; end
    if (iflag == 0); nl=nl+1; ndeg(nl)=1; end
    if (iflag == 1); ndeg(nl)=ndeg(nl)+1; end
  end
  ndeg = ndeg(1:nl);
  if(sum(ndeg)<n); error('Some problem ...'); end
  
  istop=0;
  for ib_deg = 1:nl
    istart=istop+1;
    istop=istart+ndeg(ib_deg)-1;
    Etmp1t = 0.0; Etmp2t = 0.0; Ext = 0.0; 
    for ib=istart:istop
      Etmp1t = Etmp1t + Etmp1(ib);
      Etmp2t = Etmp2t + Etmp2(ib);
      Ext = Ext + Ex(ib);
    end

    fact = ndeg(ib_deg);
    Etmp1(istart:istop) = Etmp1t / fact;
    Etmp2(istart:istop) = Etmp2t / fact;
    Ex(istart:istop) = Ext / fact;
  end

  obj.Ex = Ex;
  if ((freq_dep == 0) || (freq_dep == 1)) % Use Esex_x and Ecoh
    obj.Esex_x = Etmp1;
    obj.Ecoh   = Etmp2;
  elseif (freq_dep == 2) % Use Eres and Eint
    obj.Eres = Etmp1;
    obj.Eint = Etmp2;
  end
end % function shiftenergy