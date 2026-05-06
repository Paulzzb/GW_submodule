function rho_xalpha = get_rho_xalpha(id, param)

  
persistent firsttime N_MAX inv_rot_index R_sampling fftgrid_c 

if nargin == 1 && (ischar(id) || (isstring(id) && isscalar(id)))
  cmd = lower(string(id));
  if cmd == "reset"
    firsttime = [];
    N_MAX = [];
    inv_rot_index = [];
    R_sampling = [];
    fftgrid_c = [];
    rho_xalpha = [];
    return;
  end
end

if isempty(N_MAX)
  N_MAX = isdf.isdf_nmax();
  firsttime = true(N_MAX, 1);
  symm_data = symmetry.get();
  % We are facing the same system, so the symmetry properties are the same.
  inv_rot_index = symm_data.inv_rot_index;
  R_sampling = cell(N_MAX, 1);
  fftgrid_c = zeros(N_MAX, 3);
end

if (firsttime(id))
  firsttime(id) = false;
  isdf_data = isdf.get(id);
  R_sampling{id} = isdf_data.bundle_struct.R_grid_bundle(isdf_data.bundle_struct.sampling2bundle, :);
  fft_data = FFT.manager('get');
  fftgrid_c(id, :) = fft_data.fftgrid;
end


r_lat_data = lattice.manager('r_lat', 'get');
symm_data = symmetry.get();


iGo = param.qs(1);
iqrot = param.qs(3);

u_xalpha1 = isdf.get_u_xalpha(id, param.is, iqrot);
u_xalpha2 = isdf.get_u_xalpha(id, param.os, iqrot);
rho_xalpha = conj(u_xalpha1) .* u_xalpha2;

Go = single(r_lat_data.Ggrid_RLU(iGo, :));
if norm(Go) < 1e-6
  return;
end
inviqrot = inv_rot_index(iqrot);
invSqGo = single(Go * symm_data.rot_mtrx_RLU_G(:, :, inviqrot));
phase_shift_coarse = exp(-2 * pi * 1i * (R_sampling{id} ./ double(fftgrid_c(id, :))) * invSqGo');
% phase_shift_coarse = exp(-2 * pi * 1i * SqR_scal * Go');
% isdf_data = isdf.get(id);
% if norm(R_sampling{id} - isdf_data.R_sampling_RLU) > 1e-6
%   error('The R_sampling is not correct');
% end
rho_xalpha = rho_xalpha .* phase_shift_coarse;

end