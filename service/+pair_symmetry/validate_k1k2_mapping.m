% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19 ZZ

function report = validate_k1k2_mapping(k1k2_representation, k1k2_mapping, marked, kptbz_RLU, rot_mtrx_RLU_G, tol)
% Validate dense pair mapping against the shared-G0 criterion.
%
% This validator targets the current 3-channel mapping convention:
%   k1k2_mapping(:, :, 1) -> representative index
%   k1k2_mapping(:, :, 2) -> rotation index
%   k1k2_mapping(:, :, 3) -> G0-table index (not explicitly checked here)
% It is kept temporarily for transition/testing and should be updated when
% the new pair_symm_m.mapping(:, :, 1:4) layout becomes the only contract.

nk = size(kptbz_RLU, 1);
nsym = size(rot_mtrx_RLU_G, 3);

report = struct();
report.ok = true;
report.errors = {};
report.nk = nk;
report.nsym = nsym;
report.nrep = size(k1k2_representation, 1);
report.coverage_ratio = nnz(marked) / numel(marked);
report.max_shared_g0_error = 0.0;

if ~isequal(size(k1k2_mapping, 1), nk) || ~isequal(size(k1k2_mapping, 2), nk) || ~isequal(size(k1k2_mapping, 3), 3)
  report.ok = false;
  report.errors{end + 1} = 'k1k2_mapping has invalid shape.';
  return;
end

if any(k1k2_mapping(:, :, 1) == 0, 'all')
  report.ok = false;
  report.errors{end + 1} = 'Found unmapped representative index in k1k2_mapping(:, :, 1).';
end
if any(k1k2_mapping(:, :, 2) == 0, 'all')
  report.ok = false;
  report.errors{end + 1} = 'Found unmapped rotation index in k1k2_mapping(:, :, 2).';
end
if any(k1k2_mapping(:, :, 3) == 0, 'all')
  report.ok = false;
  report.errors{end + 1} = 'Found unmapped G0 index in k1k2_mapping(:, :, 3).';
end

for ik1 = 1:nk
  for ik2 = 1:nk
    rep_idx = int32(k1k2_mapping(ik1, ik2, 1));
    rot_idx = int32(k1k2_mapping(ik1, ik2, 2));

    if rep_idx < 1 || rep_idx > int32(size(k1k2_representation, 1))
      report.ok = false;
      report.errors{end + 1} = sprintf('Invalid representative index at (%d,%d): %d', ik1, ik2, rep_idx);
      continue;
    end

    if rot_idx < 1 || rot_idx > int32(nsym)
      report.ok = false;
      report.errors{end + 1} = sprintf('Invalid rotation index at (%d,%d): %d', ik1, ik2, rot_idx);
      continue;
    end

    rep_pair = k1k2_representation(rep_idx, :);
    k1_rep = kptbz_RLU(rep_pair(1), :);
    k2_rep = kptbz_RLU(rep_pair(2), :);
    S = rot_mtrx_RLU_G(:, :, rot_idx);

    % Reconstruct the image of the representative under S.
    k1_rot = double(k1_rep * S);
    k2_rot = double(k2_rep * S);

    target1 = double(kptbz_RLU(ik1, :));
    target2 = double(kptbz_RLU(ik2, :));

    % Accept either direct order or swapped order, because pair ordering is
    % treated as conjugate-related and not a distinct orbit.
    [ok_direct, err_direct] = local_check_shared_g0(k1_rot, k2_rot, target1, target2, tol);
    [ok_swap, err_swap] = local_check_shared_g0(k1_rot, k2_rot, target2, target1, tol);

    if ok_direct
      report.max_shared_g0_error = max(report.max_shared_g0_error, err_direct);
    elseif ok_swap
      report.max_shared_g0_error = max(report.max_shared_g0_error, err_swap);
    else
      report.ok = false;
      report.errors{end + 1} = sprintf('Shared-G0 condition failed at (%d,%d), rep=%d, rot=%d', ik1, ik2, rep_idx, rot_idx);
    end

    if k1k2_mapping(ik1, ik2, 1) ~= k1k2_mapping(ik2, ik1, 1)
      report.ok = false;
      report.errors{end + 1} = sprintf('Asymmetric representative mapping at (%d,%d)', ik1, ik2);
    end
  end
end

if ~all(marked, 'all')
  report.ok = false;
  report.errors{end + 1} = 'marked matrix is not fully covered.';
end

if report.ok
  fprintf('Validation passed: mapping is complete and shared-G0 constraints hold.\n');
else
  fprintf('Validation failed with %d issue(s).\n', numel(report.errors));
end
end

function [ok, err] = local_check_shared_g0(k1_rot, k2_rot, target1, target2, tol)
% Check whether target1,target2 are reached from k1_rot,k2_rot by one same G0.
%
% Needed condition:
%   target1 - k1_rot = G0,
%   target2 - k2_rot = G0,
% where G0 is an integer vector in reciprocal-lattice coordinates.

diff1 = target1 - k1_rot;
diff2 = target2 - k2_rot;

g1 = round(diff1);
g2 = round(diff2);

err1 = norm(diff1 - g1);
err2 = norm(diff2 - g2);
shared_err = norm(g1 - g2);

ok = (err1 <= tol) && (err2 <= tol) && (shared_err <= tol);
err = max([err1, err2, shared_err]);
end
