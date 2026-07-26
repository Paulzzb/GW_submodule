% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/08

function [wf_on_coarse, R_coarse_RLU, R_rot_coarse] = coeff_coarse_wf_filter( ...
    wf_on_coarse, R_coarse_RLU, R_rot_coarse, id)
% Drop coarse sites with tiny |wf|^2, then warn on CCH(q=0) condition number.

  amp_threshold = 1e-8;
  nmu = size(wf_on_coarse, 1);
  amp2 = sum(abs(reshape(wf_on_coarse, nmu, [])).^2, 2);

  keep_mask = amp2 >= amp_threshold;
  n_drop = nnz(~keep_mask);
  if n_drop > 0
    fprintf(['isdf.coeff.coeff_coarse_wf_filter: dropped %d/%d coarse sites ', ...
      'with ||wf||^2 < %.1e (id=%d).\n'], n_drop, nmu, amp_threshold, int32(id));
  end

  if ~any(keep_mask)
    error('isdf:coeff_coarse_wf_filter:AllDropped', ...
      'All %d coarse sites dropped (||wf||^2 < %.1e) for id=%d.', nmu, amp_threshold, int32(id));
  end

  wf_on_coarse = wf_on_coarse(keep_mask, :, :, :);
  R_coarse_RLU = R_coarse_RLU(keep_mask, :);
  R_rot_coarse = R_rot_coarse(keep_mask, :);

  local_warn_cchq0_condition(wf_on_coarse, id);
end

function local_warn_cchq0_condition(wf_on_coarse, id)
  isdf_data = isdf.get(id);
  k_data = lattice.manager('k', 'get');
  nbz = k_data.nbz;
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    return;
  end
  % warning('isdf:coeff_coarse_wf_filter:CCHq0CheckInvalid', ...
  %   'local_warn_cchq0_condition is incorrect and must not be used.');

  nrange1 = double(isdf_data.nrange1(:).');
  nrange2 = double(isdf_data.nrange2(:).');
  if isempty(nrange1) || isempty(nrange2)
    return;
  end

  if size(wf_on_coarse, 3) < 1 || size(wf_on_coarse, 4) < 1
    return;
  end

  Nisdf = size(wf_on_coarse, 1);
  Psixga = zeros(Nisdf, length(nrange1), nbz);
  Phixga = zeros(Nisdf, length(nrange2), nbz);


  
  CCHq0 = isdf.prod(Psixga, Psixga, Phixga, Phixga);
  cnd = condest(CCHq0);
  if ~isfinite(cnd) || cnd > 1e6
    warning('isdf:coeff_coarse_wf_filter:CCHq0IllConditioned', ...
      'CCHq(q=0) is ill-conditioned before adaptive update for id=%d: condest=%.6e (>1e6).', ...
      int32(id), cnd);
  end
end

