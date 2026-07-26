function vcoul = coulomb_main(GWinfo, config, iq_ibz)
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
tol_zero = 1e-7;
Godby_const = 7.44;
spherical_const = 7.7956;


% q weight
DL_vol = det(GWinfo.supercell);
RL_vol = (2*pi)^3 / DL_vol;
nbz = GWinfo.tmp_devel.nbz;
d3q_factor = RL_vol / nbz;
q_weight = d3q_factor / (2*pi)^3;

% Coulomb q=0 regularization (Godby region)
reg_q_m2 = Godby_const / (2*pi)^3 * d3q_factor^(1/3);
reg_q_m2 = reg_q_m2 * eightpi;

kpt = GWinfo.bz_samp.kpt(iq_ibz, :);
gvec = GWinfo.gvec;
qpg_RLU = gvec.components;
fftgrid = gvec.fftgrid;
nfftgridpts = gvec.nfftgridpts;

% idxnz = data.reciprocal_grid_info.idxnz{1};
% xyz = data.reciprocal_grid_info.xyz{1};
% wfncut = data.reciprocal_grid_info.wfncut;

 
% Convert to Cartesian coordinates in Ang^-1
% supercell = [a1, a2, a3];
tmp_devel = GWinfo.tmp_devel;
a1a2a3 = tmp_devel.a1a2a3;
b1b2b3 = tmp_devel.b1b2b3;


Ggrid_RLU = tmp_devel.Ggrid_RLU;
qpt_RLU = tmp_devel.kpt_RLU(iq_ibz, :);

qpg_RLU = double(qpg_RLU) + qpt_RLU;
qpg_RLU = kpt_2bz(qpg_RLU, fftgrid);
Gcart = qpg_RLU * b1b2b3';
% Gcart = double(xyz) * b1b2b3';
% Compute |q+G|^2, where q=[0 0 0] currently
if norm(qpt_RLU) < tol_zero 
  qequal0 = true;
else
  qequal0 = false;
end
qG2 = sum(Gcart.^2, 2);                 

% Calculate the coulomb potential
% ngcomb = length(Ggrid_coul.idxnz);

  switch trunc_method
    case 0  % No truncation (3D)
      vcoul = q_weight * eightpi ./ qG2;        % v(q+G) = 8π / |q+G|²
      % vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
      vcoul(qG2 < tol_zero) = reg_q_m2;      
    case 2 % 
      trunc_factor = 1-cos(sqrt(qG2) * trunc_param);
      % v(q+G) = 8π / |q+G|² * (1-cos(|q+G|*trunc_param))
      vcoul = q_weight * eightpi ./ qG2 .* trunc_factor;
      if qequal0
        vcoul(qG2 < tol_zero) = 0;           % avoid div-by-zero
      end
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