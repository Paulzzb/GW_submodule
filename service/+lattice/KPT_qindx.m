% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19 ZZ

function KPT_qindx(k, q)
  % qindx_S(ikibz, iqbz, 1) = ikpbz
  % qindx_S(ikibz, iqbz, 2) = iGo
  % kp + Go = k - q
  % 
  r_lat_m = lattice.manager('r_lat', 'get');
  fft_m = FFT.manager('get');

  nibz = k.nibz;
  nbz = q.nbz;
  
  %
  r_lat_m.qindx_S = int32( zeros(nibz, nbz, 2) );
  r_lat_m.qindx_X = int32( zeros(nibz, nbz, 2) );
  % qindx_X(iq,ikbz,1)=okbz
  % qindx_X(iq,ikbz,2)=iGo
  %
  for ikibz = 1:nibz
    kibz = k.kpt_RLU(ikibz, :);
    for iqbz = 1:nbz
      qbz = q.kptbz_RLU(iqbz, :);
      k_q = kibz - qbz;
      for ikpbz = 1:nbz
        kpbz = k.kptbz_RLU(ikpbz, :);
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
      G0 = int32( round(G0) );
      indG0 = find_gvec_in_glist( G0, r_lat_m.Ggrid_RLU, fft_m.fftgrid );
      r_lat_m.qindx_S(ikibz, iqbz, 2) = indG0;
    end
  end

  for iqibz = 1:nibz
    qibz = q.kpt_RLU(iqibz, :);
    for ikbz = 1:nbz
      kbz = k.kptbz_RLU(ikbz, :);
      k_q = kbz - qibz;
      for okbz = 1:nbz
        kbz_ok = k.kptbz_RLU(okbz, :);
        kpt_diff = k_q - kbz_ok;
        if norm( kpt_diff - round(kpt_diff) ) < r_lat_m.tol
          break;
        end
        if okbz == nbz
          error('No matching okbz found');
        end
      end
      r_lat_m.qindx_X(iqibz, ikbz, 1) = okbz;
      % 
      G0 = kbz - qibz - kbz_ok;
      G0 = int32( round(G0) );
      indG0 = find_gvec_in_glist( G0, r_lat_m.Ggrid_RLU, fft_m.fftgrid );
      r_lat_m.qindx_X(iqibz, ikbz, 2) = indG0;
    end
  end

  lattice.manager('r_lat', 'save2mod', r_lat_m);
end