function KPT_qindx(k, q)
  % qindx_S(ikibz, iqbz, 1) = ikpbz
  % qindx_S(ikibz, iqbz, 2) = iGo
  % kp + Go = k - q
  % 
  r_lat_m = lattice.manager('r_lat', 'get');

  nibz = k.nibz;
  nbz = q.nbz;
  
  %
  r_lat_m.qindx_S = int32( zeros(nibz, nbz, 2) );
  r_lat_m.qindx_X = int32( zeros(nbz, nbz, 2) );
  %
  for ikibz = 1:nibz
    kibz = k.kpt_RLU(ikibz, :);
    for iqbz = 1:nbz
      qbz = q.kptbz_RLU(iqbz, :);
      k_q = kibz - qbz;
      for ikpbz = 1:nbz
        kpbz = k.kpt_RLU(ikpbz, :);
        kpt_diff = k_q - kpbz;
        if norm( kpt_diff - round(kpt_diff) ) < r_lat_m.tol
          break;
        end
        if ikpbz == nbz
          error('No matching kpbz found');
        end
      end
      r_lat_m.qindx_S(ikibz, iqbz, 1) = ikpbz;
      % 
      G0 = kibz - qbz - kpbz;
      G0 = round(G0);
      indG0 = find_gvec_in_glist( G0, Ggrid, fftgrid );
      r_lat_m.qindx_S(ikibz, iqbz, 2) = indG0;
    end
  end

  lattice.manager('r_lat', 'save2mod', r_lat_m);
end