% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function init_from_indices(id, idx_mu)
%INIT_FROM_INDICES  Build initial ISDF bundle_struct from fine-grid indices.
%
%   isdf.rsymm.init_from_indices(id, idx_mu)
%
% Primary output: bundle_struct (fine_grid_lin, R_grid_bundle, WF_bundle,
% R_rot_in_bundle, sampling2bundle, sizes).
%
% Also dual-writes legacy point fields on isdf_m (R_sampling_RLU, coeff_seper,
% R_rot_coarse, nisdf, N_coarse, N_extra, fftgrid_c) so adaptive / SC paths
% keep working until those fields are retired.
%
% Does not set interp_scheme (set by isdf.coeff.gen_coeff).

  if nargin < 1 || isempty(id)
    output.err('ISDF id is required.');
  end
  if nargin < 2 || isempty(idx_mu)
    output.err('idx_mu must be a non-empty index vector.');
  end

  fft_data = FFT.get();
  wf_data = wave_functions.get();
  symm_data = symmetry.get();
  isdf_data = isdf.get(id);

  idx_mu = int32(idx_mu(:));
  Nmu = int32(numel(idx_mu));
  nr = double(fft_data.nr);
  if any(double(idx_mu) < 1) || any(double(idx_mu) > nr)
    output.err('idx_mu entries must lie in 1..nr (nr=%d).', nr);
  end
  if numel(unique(double(idx_mu))) ~= double(Nmu)
    output.warn('init_from_indices: idx_mu contains duplicates; keeping as-is.');
  end

  R_sampling_RLU = double(fft_data.Rgrid_RLU(idx_mu, :));

  % k-axis is IBZ (nkibz), same as wf_data.c / get_coeff / gen_bundle.
  % BZ images are obtained later via R_rot_* + WF_apply_symm, not stored here.
  nb = wf_data.nb;
  nibz = wf_data.nk;
  nspin = wf_data.nspin;
  coeff = zeros(double(Nmu), nb, nibz, nspin);
  ispin = 1;
  for ikibz = 1:nibz
    for ib = 1:nb
      coeff(:, ib, ikibz, ispin) = wf_data.c(idx_mu, ib, ikibz, ispin);
    end
  end

  % Local R_rot among selected points (identity fallback if orbit leaves the set).
  lin2loc = zeros(nr, 1, 'int32');
  lin2loc(double(idx_mu)) = int32(1:double(Nmu));
  R_rot_coarse = zeros(double(Nmu), symm_data.nsym, 'int32');
  n_open = 0;
  if isempty(fft_data.R_rot) || size(fft_data.R_rot, 2) < symm_data.nsym
    R_rot_coarse(:, :) = repmat(int32(1:double(Nmu)).', 1, symm_data.nsym);
    n_open = double(Nmu) * double(symm_data.nsym);
  else
    for isym = 1:symm_data.nsym
      rot_lin = double(fft_data.R_rot(double(idx_mu), isym));
      loc = lin2loc(rot_lin);
      miss = loc <= 0;
      n_open = n_open + nnz(miss);
      loc(miss) = int32(find(miss));
      R_rot_coarse(:, isym) = loc;
    end
  end
  if n_open > 0
    output.warn( ...
      'init_from_indices: index set is not symmetry-closed (%d mappings left the set; used identity). Prefer adaptive refine.', ...
      n_open);
  end

  % --- Primary: bundle_struct ---
  bundle_struct = struct();
  bundle_struct.N_bundle = Nmu;
  bundle_struct.N_coarse = Nmu;
  bundle_struct.N_sampling = Nmu;
  bundle_struct.R_grid_bundle = R_sampling_RLU;
  bundle_struct.R_rot_in_bundle = R_rot_coarse;
  bundle_struct.WF_bundle = coeff;
  bundle_struct.sampling2bundle = int32((1:double(Nmu)).');
  bundle_struct.fine_grid_lin = idx_mu;
  isdf_data.bundle_struct = bundle_struct;

  % --- Dual-write legacy point fields (to be retired) ---
  isdf_data.nisdf = Nmu;
  isdf_data.N_coarse = Nmu;
  isdf_data.N_extra = int32(0);
  isdf_data.coeff_seper = coeff;
  isdf_data.fftgrid_c = int32(fft_data.fftgrid(:).');
  isdf_data.R_sampling_RLU = R_sampling_RLU;
  isdf_data.R_rot_coarse = R_rot_coarse;
  isdf_data.R_rot_extra = int32(zeros(0, 0));

  isdf.save2mod(isdf_data, id);

  output.msg('rs', 'init_from_indices: id=%d  Nmu=%d', int32(id), Nmu);
end
