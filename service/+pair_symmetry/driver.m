% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/31

function driver(~, config)
  % Build pair-symmetry mapping in manager cache.
  %
  % Stored convention in pair_symm_m:
  %   representation(irep, :) = [ik1_rep, ik2_rep]
  %   mapping(ik1p, ik2p, :) encodes how a target pair (ik1p, ik2p)
  %   is generated from one representative pair under symmetry:
  %     mapping(:, :, 1) -> irep, where representation(irep, :) = [ik1_rep, ik2_rep]
  %     mapping(:, :, 2) -> isym
  %     mapping(:, :, 3) -> iG0, where g0_table(iG0, :) is the shared G0
  %
  % The intended relation is
  %   k1' = S * k1_rep + G0,
  %   k2' = S * k2_rep + G0,
  % with the same reciprocal shift G0 applied to both components.
  % WARNING: mapping convention may still evolve.

  if nargin < 2 || isempty(config)
    config = struct();
  end

  tol = 1e-5;
  if isfield(config, 'pair_symmetry_tol')
    tol = double(config.pair_symmetry_tol);
  end

  warning('pair_symmetry:UnstableModule', ...
    ['pair_symmetry.driver is initialized with an evolving mapping contract; ', ...
     'mapping construction and validation contracts may change.']);

  k_data = lattice.manager('k', 'get');
  symm_data = symmetry.manager('get');

  kptbz_RLU = double(k_data.kptbz_RLU);
  rot_mtrx_RLU_G = double(symm_data.rot_mtrx_RLU_G);
  nk = double(k_data.nbz);

  pair_data = pair_symmetry.base.pair_symm_m(nk);
  [representation, mapping, flagconj, g0_table] = ...
    build_pair_mapping(kptbz_RLU, rot_mtrx_RLU_G, tol);
  weights = build_representation_weights(mapping, size(representation, 1));

  pair_data.representation = representation;
  pair_data.nrep = int32(size(representation, 1));
  pair_data.weights = weights;
  pair_data.mapping = mapping;
  pair_data.flagconj = flagconj;
  pair_data.g0_table = g0_table;
  pair_data.assigned = true;
  pair_data.allocated = true;

  pair_symmetry.save2mod(pair_data);
end

function [representation, mapping, flagconj, g0_table] = build_pair_mapping(kptbz_RLU, rot_mtrx_RLU_G, tol)
% Enumerate representative pairs and fill the dense nk x nk mapping table.
%
% We scan the full pair grid (ik1, ik2). The first still-unassigned pair is
% taken as the next representative. Its symmetry orbit is then expanded and
% every reached target pair stores a back-reference to that representative.

nk = size(kptbz_RLU, 1);
nsym = size(rot_mtrx_RLU_G, 3);

representation = int32(zeros(0, 2));
mapping = int32(zeros(nk, nk, 3));
flagconj = false(nk, nk);
g0_table = int32(zeros(0, 3));

for ik1 = 1:nk
  for ik2 = 1:nk
    if mapping(ik1, ik2, 3) ~= 0
      continue;
    end

    % New orbit representative.
    representation(end + 1, :) = int32([ik1, ik2]); %#ok<AGROW>
    irep = int32(size(representation, 1));

    [orbit_pairs, orbit_rots, orbit_g0] = ...
      build_orbit_with_shared_g0(ik1, ik2, kptbz_RLU, rot_mtrx_RLU_G, nsym, tol);

    if isempty(orbit_pairs)
      orbit_pairs = int32([ik1, ik2]);
      orbit_rots = int32(1);
      orbit_g0 = int32([0, 0, 0]);
    end

    for ip = 1:size(orbit_pairs, 1)
      ia = double(orbit_pairs(ip, 1));
      ib = double(orbit_pairs(ip, 2));
      irot = int32(orbit_rots(ip));
      g0 = int32(orbit_g0(ip, :));

      % Store G0 by indirection so mapping remains int32-valued and compact.
      [iG0, g0_table] = find_or_add_g0(g0_table, g0);

      if mapping(ia, ib, 2) == 0
        mapping(ia, ib, 1) = irep;
        mapping(ia, ib, 2) = irot;
        mapping(ia, ib, 3) = iG0;
        flagconj(ia, ib) = false;
      end

      % Fill the swapped target as conjugate-related branch.
      % We keep the same representative/source pair but record that the
      % target ordering is reversed relative to the direct branch.
      if mapping(ib, ia, 2) == 0
        mapping(ib, ia, 1) = irep;
        mapping(ib, ia, 2) = irot;
        mapping(ib, ia, 3) = iG0;
        flagconj(ib, ia) = (ia ~= ib);
      end
    end
  end
end
end

function [orbit_pairs, orbit_rots, orbit_g0] = build_orbit_with_shared_g0(ik1, ik2, kptbz_RLU, rot_mtrx_RLU_G, nsym, tol)
% Expand one representative pair through all symmetry operations.
%
% For each S, we search BZ indices (ia, ib) such that
%   k(ia) = S*k1 + G0,
%   k(ib) = S*k2 + G0,
% with the same integer RLU shift G0 in both equations.

k1 = kptbz_RLU(ik1, :);
k2 = kptbz_RLU(ik2, :);
nk = size(kptbz_RLU, 1);

orbit_pairs = int32(zeros(0, 2));
orbit_rots = int32(zeros(0, 1));
orbit_g0 = int32(zeros(0, 3));

for irot = 1:nsym
  S = rot_mtrx_RLU_G(:, :, irot);
  k1_rot = k1 * S;
  k2_rot = k2 * S;

  [ia, ib, g0] = find_best_shared_shift_match(k1_rot, k2_rot, kptbz_RLU, nk, tol);
  if ia < 1
    continue;
  end

  orbit_pairs(end + 1, :) = int32([ia, ib]); %#ok<AGROW>
  orbit_rots(end + 1, 1) = int32(irot); %#ok<AGROW>
  orbit_g0(end + 1, :) = int32(g0); %#ok<AGROW>
end

if isempty(orbit_pairs)
  return;
end

% Multiple symmetry operations may hit the same target pair. Keep the first
% one so the mapping is double-valued.
[~, uniq_idx] = unique(double(orbit_pairs), 'rows', 'stable');
orbit_pairs = orbit_pairs(uniq_idx, :);
orbit_rots = orbit_rots(uniq_idx, :);
orbit_g0 = orbit_g0(uniq_idx, :);
end

function [best_ia, best_ib, best_g0] = find_best_shared_shift_match(k1_rot, k2_rot, kptbz_RLU, nk, tol)
% Search the BZ mesh for a pair matched by one common reciprocal shift.
%
% diff1 = k(ia) - S*k1, diff2 = k(ib) - S*k2.
% A valid match requires both diffs to be integer vectors and equal:
%   diff1 = G0, diff2 = G0.
% Among valid matches, prefer G0 = 0 and then the smallest |G0|_1.

best_ia = int32(-1);
best_ib = int32(-1);
best_g0 = int32([0, 0, 0]);
best_score = inf;

for ia = 1:nk
  diff1 = double(kptbz_RLU(ia, :) - k1_rot);
  g0_1 = round(diff1);
  if norm(diff1 - g0_1) > tol
    continue;
  end

  for ib = 1:nk
    diff2 = double(kptbz_RLU(ib, :) - k2_rot);
    g0_2 = round(diff2);
    if norm(diff2 - g0_2) > tol
      continue;
    end

    if norm(g0_1 - g0_2) > tol
      continue;
    end

    if norm(g0_1) < 0.5
      score = 0;
    else
      score = norm(g0_1, 1) + 1;
    end

    if score < best_score
      best_score = score;
      best_ia = int32(ia);
      best_ib = int32(ib);
      best_g0 = int32(g0_1);
    end
  end
end
end

function [idx, g0_table] = find_or_add_g0(g0_table, g0)
% Deduplicate shared G0 values and return the 1-based table index.
for i = 1:size(g0_table, 1)
  if isequal(g0_table(i, :), g0)
    idx = int32(i);
    return;
  end
end

g0_table(end + 1, :) = int32(g0); %#ok<AGROW>
idx = int32(size(g0_table, 1));
end

function weights = build_representation_weights(mapping, nrep)
% Compute normalized representative weights from dense pair coverage.
%
% Let irep = mapping(ik1p, ik2p, 1). The weight for representative irep is
% the fraction of target pairs mapped to it.

if nrep <= 0
  weights = double(zeros(0, 1));
  return;
end

rep_idx = double(mapping(:, :, 1));
counts = accumarray(rep_idx(:), 1, [nrep, 1]);
weights = double(counts / numel(rep_idx));
end
