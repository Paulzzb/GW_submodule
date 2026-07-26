% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function vcoul()
  coulomb_m = coulomb.get();
  r_lat_m = lattice.manager('r_lat', 'get');
  q = lattice.manager('q', 'get');

  eightpi = 8*pi;
  fourpi = 4*pi;
  Godby_const = 7.44;
  spherical_const = 7.7956;

  trunc_method = coulomb_m.trunc_method;
  trunc_param = coulomb_m.trunc_param;
  nbz = q.nbz;
  RL_vol = r_lat_m.RL_vol;
  d3q_factor = RL_vol / nbz;
  q_weight = d3q_factor / (2*pi)^3;

  reg_q_m2 = Godby_const / (2*pi)^3 * d3q_factor^(1/3);
  reg_q_m2 = reg_q_m2 * eightpi;

  coulomb_m.vcoul0 = reg_q_m2;

  for iqibz = 1:q.nibz
    % q_weight = q.weights(iqibz);
    bare_qpg = coulomb_m.bare_qpg(:, iqibz);
    switch trunc_method
      case 0
        vcoul = q_weight * eightpi ./ (bare_qpg.^2); 
      case 2
        trunc_factor = 1-cos(trunc_param*bare_qpg);
        vcoul = q_weight * eightpi .* trunc_factor ./ (bare_qpg.^2);
      otherwise
        error('Unsupported truncation method');
    end
    if iqibz == 1
      vcoul(1) = 0.0;
    end
    coulomb_m.vcoul(:, iqibz) = vcoul;
  end
  % vcoul0
  switch trunc_method
    case 0
      coulomb_m.vcoul0 = reg_q_m2;
    case 2
      coulomb_m.vcoul0 = q_weight * fourpi * trunc_param.^2;
    otherwise
      error('Unsupported truncation method');
  end
  %
  coulomb.save2mod(coulomb_m);
end
