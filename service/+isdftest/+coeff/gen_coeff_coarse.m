% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function gen_coeff_coarse(id)
% ISDF coarse index / grid API (dispatcher, no nested functions).
%
%   [Nmu, ind_mu] = isdftest.coeff.gen_coeff_coarse()
%       R-grid subset selection (gcd / step); Nmu_tmp = isdf_ratio * nrep.
%
%   [fftgrid_i, fft_sz, nfft_d, fftgrid_c, Nmu, R_coarse_RLU, R_rot_coarse] = ...
%       isdftest.coeff.gen_coeff_coarse('fft_grid')
%       Regular coarse FFT box; Nmu_tmp = isdf_ratio * nb.
%
%   wf_on_coarse = isdftest.coeff.gen_coeff_coarse('wf', fftgrid_i, fft_sz, fftgrid_c, Nmu, R_coarse_RLU)
%       Map wavefunctions to coarse RLU sites.

  fft_data = FFT.get();

  isdftest.coeff.coeff_coarse_fft_grid_build(id);
  isdf_data = isdftest.get(id);
  fftgrid_c = isdf_data.fftgrid_c;
  R_coarse_RLU = isdf_data.R_sampling_RLU;
  R_rot_coarse = isdf_data.R_rot_coarse;

  wf_on_coarse = isdftest.coeff.coeff_coarse_wf_extract(fftgrid_c, R_coarse_RLU, id);

  % Rewrite isdf_data (nisdf must match coeff rows; same as prod(fftgrid_c) from build)
  isdf_data.nisdf = int32(size(wf_on_coarse, 1));
  isdf_data.N_coarse = int32(size(R_rot_coarse, 1));
  isdf_data.N_extra = int32(0);
  isdf_data.coeff_seper = wf_on_coarse;
  scal =  double(fft_data.fftgrid) ./ double(fftgrid_c);
  isdf_data.R_sampling_RLU = R_coarse_RLU .* scal;
  isdf_data.interp_scheme = "coarse";
  isdf_data.R_rot_coarse = R_rot_coarse;
  %
  bundle_struct = struct();
  bundle_struct.N_bundle = isdf_data.N_coarse;
  bundle_struct.N_coarse = isdf_data.N_coarse;
  bundle_struct.N_sampling = isdf_data.N_coarse;
  bundle_struct.R_grid_bundle = R_coarse_RLU .* scal;
  bundle_struct.R_rot_in_bundle = R_rot_coarse;
  bundle_struct.WF_bundle = wf_on_coarse;
  bundle_struct.sampling2bundle = 1:isdf_data.N_coarse;
  isdf_data.bundle_struct = bundle_struct;
  %
  isdftest.save2mod(isdf_data, id);
  
end
