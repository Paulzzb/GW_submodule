function vcoul_cell = construct_vcoul_k(data, config, gvec_list)
% construct_vcoul - construct the electrostatic potential from the charge density
% config.CUTOFFS.coulomb_truncation_method = 
%   = 0,  no truncation (3D)
%   = 2,  0D spherical truncation
%   = 4,  1D cell wire truncation
%   = 5,  0D cell box truncation
%   = 6,  2D slab truncation
%   = 7,  supercell truncation (3D), experimental

eightpi = 8*pi;
fourpi = 4*pi;
tol_zero = 1e-10;

trunc_method = config.CUTOFFS.coulomb_truncation_method;
trunc_param = config.CUTOFFS.coulomb_truncation_parameter;

supercell = data.sys.supercell;
recip_lattice = 2*pi*inv(supercell'); % rows are b1, b2, b3

nkibz = data.nkibz;
qpoint_list = data.kibz;
vcoul_cell = cell(nkibz, 1);

for ik = 1:nkibz
  gvec = gvec_list;
  qpoint = qpoint_list(ik, :);
  Gcart = gvec.components * recip_lattice;
  % Compute |q+G|^2, where q=[0 0 0] currently
  qG = qpoint + Gcart;
  qG2 = sum(qG.^2, 2);                 
   
  % Calculate the coulomb potential

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
  vcoul_cell{ik} = vcoul;
end % for ik


end % EOF