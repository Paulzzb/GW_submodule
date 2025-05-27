function qpenergy = getEqp(qpenergy)
  
  if ((qpenergy.freq_dep == 0) || (qpenergy.freq_dep == 1)) % Use Esex_x and Ecoh
    Etmp1 = qpenergy.Esex_x;
    Etmp2 = qpenergy.Ecoh;
  elseif (qpenergy.freq_dep == 2) % Use Eres and Eint
    Etmp1 = qpenergy.Eres;
    Etmp2 = qpenergy.Eint;
  end

  qpenergy.Sig = qpenergy.Ex + Etmp1 + Etmp2;
  qpenergy.Eqp = qpenergy.Sig - qpenergy.Vxc + qpenergy.ev;
end