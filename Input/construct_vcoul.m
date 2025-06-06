% function vcoul = construct_vcoul(mill, supercell, amin, truncation)
function vcoul = construct_vcoul(data, config)

trun_method = config.CUTOFFS.coulomb_truncation_method;
cutoff = config.CUTOFFS.coulomb_cutoff;
eightpi = 8*pi;
fourpi = 4*pi;
tol_zero = 1e-10;


idxnz = data.reciprocal_grid_info.idxnz;
xyz = data.reciprocal_grid_info.xyz;
wfncut = data.reciprocal_grid_info.wfncut;
 
% Convert to Cartesian coordinates in Ang^-1
supercell = data.sys.supercell;
recip_lattice = 2*pi*inv(supercell'); % rows are b1, b2, b3
Gcart = xyz * recip_lattice;
% Compute |q+G|^2, where q=[0 0 0] currently
qG2 = sum(Gcart.^2, 2);                 
 
% Select desired reciprocal vectors under coulomb_cutoff if coulomb_cutoff < wfncut
if (cutoff < wfncut)
  Ggrid_coul= [];
  % Filter where G² <= cutoff
  new_idx = find(G2 <= cutoff);
  Ggrid_coul.xyz = xyz(new_idx,:);
  Ggrid_coul.idxnz = idxnz(new_idx,:);
  Gcart = Gcart(new_idx,:);
  qG2 = qG2(new_idx,:);
else
  Ggrid_coul = data.reciprocal_grid_info; 
end

% Calculate the coulomb potential
ngcomb = length(Ggrid_coul.idxnz);

switch trun_method
  case 0  % No truncation (3D)
    vcoul = eightpi ./ qG2;        % v(q+G) = 8π / |q+G|²
    vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero

  case 4  % Wire truncation (1D)
    % Apply wire-truncated Coulomb potential
    % Formula (e.g., from Rozzi et al. 2006):
    % v(q+G) = 4π / |q+G|² * [1 - exp(-|q+G|*L) * (|q+G|*L + 1)]
    L = truncval(1); % length along non-periodic directions
    qGnorm = sqrt(qG2);
    trunc_factor = 1 - exp(-qGnorm * L) .* (qGnorm * L + 1);
    vcoul = fourpi ./ qG2 .* trunc_factor;
    vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero

  case 5  % Box truncation (0D)
    % Use erfc-based cutoff (example form)
    L = truncval(1);
    qGnorm = sqrt(qG2);
    vcoul = fourpi ./ qG2 .* (1 - exp(-qGnorm.^2 * L^2));
    vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero

  case 6  % Slab truncation (3D)
    % Slab truncation along z (say Gz = G(:,3))
    Lz = truncval(3);
    Gz = qG(:,3);
    vcoul = eightpi ./ qG2 .* (1 - exp(-abs(Gz) * Lz));
    vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero

  otherwise
    error('Unsupported truncation type: %d', trunc_type);
end
end;



end % EOF