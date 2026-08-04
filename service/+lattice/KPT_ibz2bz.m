% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function k = KPT_ibz2bz(k)
% KPT_ibz2bz - expand irreducible k-points in the first Brillouin zone
% to the full k-points set (BZ)
% Then construct kptbz
  
  symm_m = symmetry.manager('get');
  r_lat_m = lattice.manager('r_lat', 'get');
  b1b2b3 = r_lat_m.b1b2b3;

  k.kptbz_RLU = double( zeros(k.nbz, 3) );
  k.kptbz_Cart = double( zeros(k.nbz, 3) );
  for ikbz = 1:k.nbz
    ikibz = k.bz2ibz(ikbz);
    isym = k.bz2rot(ikbz);
    rot_Cart = double( symm_m.rot_mtrx_Cart(:, :, isym) );
    kibz_Cart = double( k.kpt_Cart(ikibz, :) );
    k.kptbz_Cart(ikbz, :) = kibz_Cart * rot_Cart;
  end
  k.kptbz_RLU = k.kptbz_Cart / b1b2b3';
  
end
