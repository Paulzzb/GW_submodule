% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/31

classdef pair_symm_m
  properties
    % Standard service-module lifecycle flags.
    assigned(1, 1) logical = false
    allocated(1, 1) logical = false

    % Number of full-BZ k-points used to define the pair grid.
    nk(1, 1) {mustBeInteger, mustBeNonnegative} = int32(0)

    % representation(irep, :) = [ik1_rep, ik2_rep].
    % Each row is one chosen source pair whose symmetry orbit covers part
    % of the full nk x nk target-pair table.
    representation(:, 2) int32 = int32(zeros(0, 2))
    nrep (1, 1) int32 = int32(0)

    % weights(irep) is the normalized weight of representative irep on the
    % full pair grid. It is typically orbit_size(irep) / (nk * nk).
    weights(:, 1) double = double(zeros(0, 1))

    % mapping(ik1p, ik2p, :):
    %   (:, :, 1) -> representative index irep, with
    %                representation(irep, :) = [ik1_rep, ik2_rep]
    %   (:, :, 2) -> rotation index (isym)
    %   (:, :, 3) -> shared G0 index
    % For a target pair (ik1p, ik2p), reconstruction uses
    %   [ik1_rep, ik2_rep] = representation(mapping(ik1p, ik2p, 1), :),
    %   k1p = S*k1_rep + G0,
    %   k2p = S*k2_rep + G0,
    % with S given by mapping(:, :, 2) and G0 by g0_table(mapping(:, :, 3), :). 
    mapping(:, :, :) int32 = int32(zeros(0, 0, 3))

    % g0_table(iG0, :) stores integer RLU triplet [g1, g2, g3]
    % corresponding to mapping(:, :, 3) == iG0.
    g0_table(:, 3) int32 = int32(zeros(0, 3))

    % flagconj(ik1p, ik2p) = true if this target pair is filled through the
    % swapped-order branch relative to the direct mapped pair.
    flagconj(:, :) logical = false(0, 0)
  end

  methods
    function obj = pair_symm_m(nk)
      if nargin == 0
        return
      end

      validateattributes(nk, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'nk');

  % Allocate dense tables for the full target pair grid. Representation
  % and g0_table are filled later by the driver.
      obj.nk = int32(nk);
      obj.representation = int32(zeros(0, 2));
      obj.weights = double(zeros(0, 1));
      obj.mapping = int32(zeros(nk, nk, 3));
      obj.g0_table = int32(zeros(0, 3));
      obj.flagconj = false(nk, nk);
      obj.allocated = true;
    end
  end
end
