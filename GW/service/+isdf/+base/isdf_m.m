% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/30

classdef isdf_m
  properties
    assigned(1, 1) logical = false
    allocated(1, 1) logical = false

    id(1, 1) {mustBeInteger} = int32(0)
    desc(1, 1) string = ""

    dbroot(1, 1) string = ""
    sys(1, 1) string = ""

    interp_scheme(1, 1) string = "" % describes the interpolation scheme to use
    nisdf(1, 1) {mustBeInteger} = int32(0)
    R_sampling_RLU(:, :) single = single(zeros(0, 3))
    % Stored as double arrays; values may be complex (avoid classdef 'complex' = zeros(...) on older MATLAB).
    coeff_seper(:, :, :, :) = zeros(0, 0, 0, 0) % nisdf * nb * nkibz * nspin
    tildeVq(:, :, :, :) = zeros(0, 0, 0, 0) % nisdf * nisdf * nkibz * nspin
    helperqG(:, :, :, :) = zeros(0, 0, 0, 0) % ng * nisdf * nkibz * nspin
    % Case 1: interp_scheme = 'coarse', Coarse FFT box
    % Case 2: interp_scheme = 'adaptive', first coarse, then a subset of fine grid points
    fftgrid_c(:, :) {mustBeInteger} = int32(zeros(0, 0))
      % for getting wavefunctions
    R_rot_coarse(:, :) {mustBeInteger} = int32(zeros(0, 0))
      % for applying symmetry to coarse ( and possibly extra) grid points
    R_rot_extra(:, :) {mustBeInteger} = int32(zeros(0, 0))
    N_coarse(1, 1) {mustBeInteger} = int32(0)
    N_extra(1, 1) {mustBeInteger} = int32(0) % N_extra + N_coarse = nisdf
    tmp
    bundle_struct
  end

  methods
    function obj = isdf_m()
      obj.allocated = true;
      obj.assigned = true;
    end
  end
end
