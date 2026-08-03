% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/13 ZZ

function isdf_exclude_point(id)
%Because we realize that, the given fftgrid could lead to ill-conditioned
% Gram--matrix, indicating removing some ill-scaled points could be
% helpful.
  rel_tol = double(1e-5);

  d_lat_data = lattice.manager('d_lat', 'get');
  DL_vol = d_lat_data.DL_vol;
  fft_data = FFT.get();
  nr = fft_data.nr;

  abs_tol = double(nr) / sqrt(DL_vol) * rel_tol;

  isdf_data = isdf.get(id);
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('isdf_exclude_point:nrange', ...
      'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', ...
      int32(id));
  end
  nrange1 = double(isdf_data.nrange1);
  nrange2 = double(isdf_data.nrange2);
  Psixga = isdf_data.coeff_seper(:, nrange1, 1, 1);
  Phixga = isdf_data.coeff_seper(:, nrange2, 1, 1);

  Nisdf_tmp = size(Psixga, 1);

  w = zeros(Nisdf_tmp, 1);
  for i = 1:Nisdf_tmp
    w(i) = isdf.adaptive.isdf_prod(Psixga(i, :), Psixga(i, :), Phixga(i, :), Phixga(i, :));
  end

  well_ind = find(w > abs_tol);
  ill_ind = find(w <= abs_tol);
  if isempty(ill_ind)
    fprintf('isdf_exclude_point: id=%d 鈥?no centroids removed (all w > abs_tol=%.6e).\n', ...
      id, double(abs_tol));
    return
  end
  if isempty(well_ind)
    error('isdf_exclude_point:AllRemoved', ...
      'All %d ISDF centroids fail w > abs_tol (abs_tol=%g); check grid/tolerance.', ...
      Nisdf_tmp, double(abs_tol));
  end

  % --- simple console report: removed centroid rows vs weight w ---
  fprintf('\n=== isdf_exclude_point (ISDF id=%d, desc=%s) ===\n', id, char(isdf_data.desc));
  fprintf('  threshold: abs_tol = %.6e (rel_tol=%.3e, nr=%.6g, sqrt(DL_vol)=%.6e)\n', ...
    double(abs_tol), double(rel_tol), double(nr), sqrt(double(DL_vol)));
  fprintf('  centroids: total %d, kept %d, removed %d\n', ...
    Nisdf_tmp, numel(well_ind), numel(ill_ind));
  fprintf('  removed row (1-based in coeff_seper) | w (diag proxy)\n');
  for k = 1:numel(ill_ind)
    ir = ill_ind(k);
    fprintf('  %8d  |  %.6e\n', ir, w(ir));
  end
  fprintf('=== end isdf_exclude_point ===\n\n');

  Nnew = numel(well_ind);
  isdf_data.nisdf = int32(Nnew);
  isdf_data.coeff_seper = isdf_data.coeff_seper(well_ind, :, :, :);
  isdf_data.R_rot_coarse = isdf_data.R_rot_coarse(well_ind, :);
  isdf_data.R_sampling_RLU = isdf_data.R_sampling_RLU(well_ind, :);
  % isdf_data.fftgrid_c = isdf_data.fftgrid_c(well_ind, :);

  % Row-aligned with centroids (same row count as coeff_seper before filter = Nisdf_tmp).
  if size(isdf_data.R_sampling_RLU, 1) == Nisdf_tmp
    isdf_data.R_sampling_RLU = isdf_data.R_sampling_RLU(well_ind, :);
  end

  if ~isempty(isdf_data.tildeVq) && size(isdf_data.tildeVq, 1) == Nisdf_tmp ...
      && size(isdf_data.tildeVq, 2) == Nisdf_tmp
    isdf_data.tildeVq = isdf_data.tildeVq(well_ind, well_ind, :, :);
  end

  if ~isempty(isdf_data.helperqG) && size(isdf_data.helperqG, 2) == Nisdf_tmp
    isdf_data.helperqG = isdf_data.helperqG(:, well_ind, :, :);
  end

  % save2mod(data, id) 鈥?not (id, data); do not assign to package name isdf.*
  
  isdf.save2mod(isdf_data, id);
end
