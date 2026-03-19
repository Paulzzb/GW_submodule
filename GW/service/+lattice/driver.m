% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function driver(data, config)
  
  symm_m = symmetry.manager('get');
  nsym = symm_m.nsym;

  % Initial setup for d_lat and r_lat
  lattice.manager('r_lat', 'free');
  lattice.manager('d_lat', 'free');
  d_lat_m = lattice.d_lattice_m();
  r_lat_m = lattice.r_lattice_m();
  d_lat_m.a1a2a3 = data.sys.supercell';
  d_lat_m.DL_vol = det(d_lat_m.a1a2a3);
  b1b2b3 = 2*pi*inv(r_lat_m.a1a2a3');
  r_lat_m.b1b2b3 = b1b2b3;
  r_lat_m.RL_vol = (2*pi)^3 / d_lat_m.DL_vol;
  lattice.manager('r_lat', 'save2mod', r_lat_m);
  lattice.manager('d_lat', 'save2mod', d_lat_m);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Lets construct G-grid
  % 
  ecut = config.CUTOFF.coulomb_cutoff;
  lattice.RL_shell_construct(ecut); 
  % 
  % Construct G_rot
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct kpt first
  kibz_RLU = data.kibz;
  nibz = size(kibz_RLU, 1);
  kweight = data.kweight;

  k = bz_samp(nibz, nsym);

  k.kweight = kweight;
  k.kpt_RLU = kibz_RLU; 
  k.kpt_Cart = k.kpt_RLU * b1b2b3';

  % Set all properties in k with respect to 'full bz'
  k = KPT_expand(k);
  q = k;
  %
  % Construct qindx*
  % 
  KPT_qindx(k, q);
  %
  lattice.manager('k', 'save2mod', k);
  lattice.manager('q', 'save2mod', q);
end