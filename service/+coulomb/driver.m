% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function driver(data, config)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Load
  symm_m = symmetry.manager('get');
  nsym = symm_m.nsym;

  % Initial setup for d_lat and r_lat
  r_lat = lattice.manager('r_lat', 'get');
  q = lattice.manager('q', 'get');
  
  ng = r_lat.ng;
  nqibz = q.nibz;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initial setup for coulomb
  coulomb_m = coulomb.base.coulomb_m(ng, nqibz);
  %
  coulomb_m.trunc_method = config.CUTOFFS.coulomb_truncation_method;
  coulomb_m.trunc_param = config.CUTOFFS.coulomb_truncation_parameter;
  %
  % Construct bare_qpg
  Ggrid_Cart = r_lat.Ggrid_Cart;
  %
  for iqibz = 1:nqibz
    qpt = q.kpt_Cart(iqibz, :);
    qgp_Cart = Ggrid_Cart + qpt;
    qgpabs = sqrt( sum(qgp_Cart.^2, 2) );
    coulomb_m.bare_qpg(:, iqibz) = qgpabs;
  end
  coulomb.save2mod(coulomb_m);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct vcoul and vcoul0 
  coulomb.vcoul();
end 