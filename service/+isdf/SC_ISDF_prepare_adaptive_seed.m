% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11

function SC_ISDF_prepare_adaptive_seed(id)
%SC_ISDF_PREPARE_ADAPTIVE_SEED  Prepare SC-replicated ISDF as adaptive initial set.
%
% Fills coeff on the supercell fine grid and sets coarse-grid metadata
% (R_rot_coarse, fine_grid_lin, interp_scheme='coarse') for adaptiveisdf.

  isdf.coeff.gen_coeff_from_fine_grid(id);
  isdf_data = isdf.get(id);
  Nmu = double(isdf_data.nisdf);
  if Nmu < 1
    error('isdf:SC_ISDF_prepare_adaptive_seed:Empty', ...
      'ISDF id=%d has nisdf < 1 after SC_ISDF.', int32(id));
  end

  if isfield(isdf_data, 'bundle_struct') && isstruct(isdf_data.bundle_struct) ...
      && isfield(isdf_data.bundle_struct, 'fine_grid_lin') ...
      && ~isempty(isdf_data.bundle_struct.fine_grid_lin)
    lin = int32(isdf_data.bundle_struct.fine_grid_lin(:));
  else
    lin = isdf.SC_ISDF_r_sampling_to_lin(isdf_data);
  end

  [lin, ia] = unique(lin, 'stable');
  if numel(ia) < Nmu
    fprintf('[SC_ISDF] unique fine_grid_lin: %d -> %d\n', Nmu, numel(ia));
    isdf_data = sc_subset_isdf_rows(isdf_data, ia, lin);
    Nmu = numel(ia);
  end

  isdf_data.N_coarse = int32(Nmu);
  isdf_data.N_extra = int32(0);
  isdf_data.R_rot_coarse = sc_build_R_rot_from_lin(lin);
  isdf_data.R_rot_extra = isdf_data.R_rot_coarse;
  isdf_data.interp_scheme = "coarse";
  isdf_data.tmp = isdf_data.coeff_seper;

  bs = isdf_data.bundle_struct;
  if ~isstruct(bs)
    bs = struct();
  end
  bs.N_coarse = int32(Nmu);
  bs.N_sampling = int32(Nmu);
  bs.N_bundle = int32(Nmu);
  bs.fine_grid_lin = lin;
  bs.sampling2bundle = int32((1:Nmu).');
  bs.R_grid_bundle = isdf_data.R_sampling_RLU(1:Nmu, :);
  bs.R_rot_in_bundle = isdf_data.R_rot_coarse;
  bs.WF_bundle = isdf_data.coeff_seper;
  isdf_data.bundle_struct = bs;

  isdf.save2mod(isdf_data, id);
  isdf.coeff.isdf_apply_symm_on_coarse('reset');
  fprintf('[SC_ISDF] prepared adaptive seed on id=%d, nisdf=%d, desc=%s\n', ...
    int32(id), int32(Nmu), char(string(isdf_data.desc)));
end

function R_rot = sc_build_R_rot_from_lin(lin_sc)
  lin_sc = int32(lin_sc(:));
  Ntot = numel(lin_sc);
  symm_data = symmetry.get();
  nsym = double(symm_data.nsym);
  fft_data = FFT.get();
  fftgrid = double(fft_data.fftgrid(:)).';
  Rgrid = double(fft_data.Rgrid_RLU);
  R_rot = zeros(Ntot, nsym, 'int32');
  for is = 1:nsym
    mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
    M2_r_RLU = Rgrid(lin_sc, :) * mtrx_RLU_R;
    if max(abs(M2_r_RLU - round(M2_r_RLU)), [], 'all') > 1e-3
      error('isdf:SC_ISDF_prepare_adaptive_seed:NonIntegerRot', ...
        'Non-integer rotation map on supercell FFT grid (symmetry %d).', is);
    end
    iv_mod = int32(mod(round(M2_r_RLU) + fftgrid, fftgrid));
    i4 = 1 + iv_mod(:, 1) + iv_mod(:, 2) * fftgrid(1) ...
      + iv_mod(:, 3) * fftgrid(1) * fftgrid(2);
    R_rot(:, is) = int32(i4);
  end
end

function isdf_data = sc_subset_isdf_rows(isdf_data, keep, lin)
  isdf_data.nisdf = int32(numel(keep));
  isdf_data.coeff_seper = isdf_data.coeff_seper(keep, :, :, :);
  isdf_data.R_sampling_RLU = isdf_data.R_sampling_RLU(keep, :);
  if ~isempty(isdf_data.R_rot_extra)
    isdf_data.R_rot_extra = isdf_data.R_rot_extra(keep, :);
  end
  if ~isempty(isdf_data.R_rot_coarse)
    isdf_data.R_rot_coarse = isdf_data.R_rot_coarse(keep, :);
  end
  isdf_data.N_coarse = int32(numel(keep));
  isdf_data.N_extra = int32(0);
  if ~isfield(isdf_data, 'bundle_struct') || ~isstruct(isdf_data.bundle_struct)
    isdf_data.bundle_struct = struct();
  end
  isdf_data.bundle_struct.fine_grid_lin = int32(lin(:));
end
