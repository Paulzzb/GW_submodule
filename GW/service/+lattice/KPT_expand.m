% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function k = KPT_expand(k)
  % KPT_expand(k)
  %   For each irreducible k-point, find all symmetrically equivalent
  %   k-points in the full BZ.

  symm_m = symmetry.manager('get');
  nsym = symm_m.nsym;
  r_lat_m = lattice.manager('r_lat', 'get');
  b1b2b3 = r_lat_m.b1b2b3;

  nibz = k.nibz;
  k.nstar = int32( zeros(nibz, 1) );
  k.star  = int32( zeros(nibz, nsym) );
  kstar = zeros(nsym, 3);
  for ikibz = 1:nibz
    kibz_Cart = k.kpt_Cart(ikibz, :);
    for isym = 1:nsym
      rot_Cart = symm_m.rot_mtrx_Cart(:, :, isym);
      Skibz = single( kibz_Cart * rot_Cart );
      Skibz = bz2ibz( 'Cart', Skibz, b1b2b3);
      Skibz_RLU = Skibz / b1b2b3';
      kstar(isym, :) = Skibz_RLU;
      k_found = false;
      for i_star = 1: k.nstar(ikibz)
        tmp = kstar(isym, :) - kstar(k.star(ikibz, i_star), :);
        diff = tmp - round(tmp);
        if norm( diff ) < 1e-5 * r_lat_m.RL_vol
          k_found = true;
          break;
        end
      end
      if ~k_found
        k.nstar(ikibz) = k.nstar(ikibz) + 1;
        k.star(ikibz, k.nstar(ikibz)) = isym;
      end
    end
  end
  
  k.nbz = sum(k.nstar);
  weights = single( k.nstar(:) ) / single( k.nbz );

  k.sstar = int32( zeros(k.nbz, 2) );
  k.kptbz_RLU   = single(zeros(k.nbz, 3));
  k.kptbz_Cart = single(zeros(k.nbz, 3));
  k.bz2ibz = int32(zeros(k.nbz, 1));
  k.bz2rot = int32(zeros(k.nbz, 1));
  
  if norm(weights - k.weights) > 1e-5
    warning('The calculated weights differ from the input weights. Using the calculated weights.');
  end
  k.weights = weights;

  k.bz2ibz = int32( zeros(k.nbz, 1) );
  k.bz2rot = int32( zeros(k.nbz, 1) );
  
  ikbz = int32( 0 );
  for ikibz = 1:nibz
    for istar = 1:k.nstar(ikibz)
      ikbz = ikbz+1;
      if istar == 1
        k.ibz2bz( ikibz ) = ikbz;
      end
      k.bz2ibz( ikbz ) = ikibz;
      k.bz2rot( ikbz ) = k.star(ikibz, istar);
    end
  end

  k = lattice.KPT_ibz2bz(k); 


end