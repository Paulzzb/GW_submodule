function vcoul = construct_vcoul(data, config, gvec)
% construct_vcoul - construct the electrostatic potential from the charge density
% config.CUTOFFS.coulomb_truncation_method = 
%   = 0,  no truncation (3D)
%   = 2,  0D spherical truncation
%   = 4,  1D cell wire truncation
%   = 5,  0D cell box truncation
%   = 6,  2D slab truncation
%   = 7,  supercell truncation (3D), experimental
% 
trunc_method = config.CUTOFFS.coulomb_truncation_method;
trunc_param = config.CUTOFFS.coulomb_truncation_parameter;
eightpi = 8*pi;
fourpi = 4*pi;
tol_zero = 1e-10;


idxnz = data.reciprocal_grid_info.idxnz{1};
xyz = data.reciprocal_grid_info.xyz{1};
wfncut = data.reciprocal_grid_info.wfncut;
 
% Convert to Cartesian coordinates in Ang^-1
supercell = data.sys.supercell;
recip_lattice = 2*pi*inv(supercell'); % rows are b1, b2, b3
Gcart = double(xyz) * recip_lattice;
Gcart = double(xyz) * recip_lattice;
% Compute |q+G|^2, where q=[0 0 0] currently
qG2 = sum(Gcart.^2, 2);                 

% Calculate the coulomb potential
% ngcomb = length(Ggrid_coul.idxnz);

  switch trunc_method
    case 0  % No truncation (3D)
      vcoul = eightpi ./ qG2;        % v(q+G) = 8π / |q+G|²
      vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
    case 2 % 
      trunc_factor = 1-cos(sqrt(qG2) * trunc_param);
      vcoul = eightpi ./ qG2 .* trunc_factor;        % v(q+G) = 8π / |q+G|² * (1-)
      vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
    case 4  % Wire truncation (1D)
      % Apply wire-truncated Coulomb potential
      % Formula (e.g., from Rozzi et al. 2006):
      % v(q+G) = 4π / |q+G|² * [1 - exp(-|q+G|*L) * (|q+G|*L + 1)]
      msg = sprintf(['GW:construct_vcoul:trunc_method', 'Wire truncation not implemented yet']);
      QPerror(msg);
      % error('GW:construct_vcoul:trunc_method', 'Wire truncation not implemented yet');
      L = truncval(1); % length along non-periodic directions
      qGnorm = sqrt(qG2);
      trunc_factor = 1 - exp(-qGnorm * L) .* (qGnorm * L + 1);
      vcoul = fourpi ./ qG2 .* trunc_factor;
      vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
  
    case 5  % Box truncation (0D)
      % Use erfc-based cutoff (example form)
      msg = sprintf(['GW:construct_vcoul:trunc_method', 'Box truncation not implemented yet']);
      QPerror(msg);
  
      error('GW:construct_vcoul:trunc_method', 'Box truncation not implemented yet');
      L = truncval(1);
      qGnorm = sqrt(qG2);
      vcoul = fourpi ./ qG2 .* (1 - exp(-qGnorm.^2 * L^2));
      vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
  
    case 6  % Slab truncation (3D)
      % Slab truncation along z (say Gz = G(:,3))
      msg = sprintf(['GW:construct_vcoul:trunc_method', 'Slab truncation not implemented yet']);
      QPerror(msg);
      error('GW:construct_vcoul:trunc_method', 'Slab truncation not implemented yet');
      Lz = truncval(3);
      Gz = qG(:,3);
      vcoul = eightpi ./ qG2 .* (1 - exp(-abs(Gz) * Lz));
      vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
    case 7
      msg = sprintf(['GW:construct_vcoul:trunc_method', 'Supercell truncation not implemented yet']);
      QPerror(msg);
      error('GW:construct_vcoul:trunc_method', 'Supercell truncation not implemented yet');
  
    otherwise
      msg = sprintf(['Unsupported truncation type: %d', trunc_method]);
      QPerror(msg);
      error('Unsupported truncation type: %d', trunc_method);
  end



end % EOF