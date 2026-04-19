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
  d_lat_m = lattice.base.d_lattice_m();
  r_lat_m = lattice.base.r_lattice_m();
  d_lat_m.a1a2a3 = data.sys.supercell';
  d_lat_m.DL_vol = det(d_lat_m.a1a2a3);
  b1b2b3 = 2*pi*inv(d_lat_m.a1a2a3');
  r_lat_m.b1b2b3 = b1b2b3;
  r_lat_m.RL_vol = (2*pi)^3 / d_lat_m.DL_vol;

  xyz = [];
  if isfield(data, 'xyz') && ~isempty(data.xyz)
    xyz = data.xyz;
  elseif isfield(data, 'sys') && isfield(data.sys, 'xyzlist') && ~isempty(data.sys.xyzlist)
    xyz = data.sys.xyzlist;
  end
  if ~isempty(xyz)
    xyz = single(xyz);
    nat = size(xyz, 1);
    d_lat_m.atom_pos = reshape(xyz, nat, 1, 3);
    if isfield(data, 'atom_symbol') && numel(data.atom_symbol) == nat
      d_lat_m.atom_symbol = data.atom_symbol(:);
    else
      d_lat_m.atom_symbol = repmat({''}, nat, 1);
    end
  end

  lattice.manager('r_lat', 'save2mod', r_lat_m);
  lattice.manager('d_lat', 'save2mod', d_lat_m);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Lets construct G-grid
  % 
  ecut = config.CUTOFFS.coulomb_cutoff;
  lattice.RL_shell_construct(ecut); 
  % 
  % Construct G_rot
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct kpt first
  kibz_Cart = data.kibz;
  nibz = size(kibz_Cart, 1);
  kweights = data.kweight;

  k = lattice.base.bz_samp_m(nibz, nsym);

  k.weights = kweights;
  k.kpt_Cart = single(kibz_Cart);
  k.kpt_RLU = k.kpt_Cart / b1b2b3';

  % Set all properties in k with respect to 'full bz'
  k = lattice.KPT_expand(k);
  q = k;
  %
  % Construct qindx*
  % 
  lattice.KPT_qindx(k, q);
  FFT.FFT_G_table();
  %
  lattice.manager('k', 'save2mod', k);
  lattice.manager('q', 'save2mod', q);
end