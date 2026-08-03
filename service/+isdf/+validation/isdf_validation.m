% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [Esum2, EsumISDF2, DiffEsum2] = isdf_validation(id, outDir)
% Run coarse ISDF energy validation for pool id (same physics as isdftest/isdf_coarse_validate_energies).
%
% Requires: isdf.coeff.gen_coeff_coarse(id) and isdf.gen_tildeVq(id) so that
%   interp_scheme == "coarse", coeff_seper holds wf_on_coarse, tildeVq is filled.
%
% Optional id: omit to use isdf.isdf_current().
% Optional outDir: directory for isdf_validate_HF_id<id>.txt (default pwd).

  if nargin < 1 || isempty(id)
    id = isdf.isdf_current();
  end
  if nargin < 2
    outDir = '';
  end

  isdf_data = isdf.get(id);

  switch isdf_data.interp_scheme
    case "coarse"
      interp_scheme = "coarse";
    case {"adaptive", "adaptive_sc", "coarse_sc", "supercell"}
      interp_scheme = "adaptive";
    otherwise
      error('isdf:isdf_validation:InvalidInterpScheme', ...
        'Invalid interp_scheme: %s', isdf_data.interp_scheme);
  end

  if isempty(isdf_data.tildeVq)
    error('isdf:isdf_validation:NoTildeVq', ...
      'tildeVq is empty; run isdf.gen_tildeVq(id) first. id=%d', id);
  end

  Nisdf = isdf_data.nisdf;
  % tildeVq = isdf_data.tildeVq;
  % if ndims(tildeVq) > 3
  %   tildeVq = tildeVq(:, :, :, 1);
  % end
  % tildeVq = squeeze(tildeVq);
  % if ndims(tildeVq) ~= 3
  %   error('isdf:isdf_validation:BadTildeVq', 'tildeVq must be 3-D (Nisdf x Nisdf x nibz) after squeeze.');
  % end

  wf_on_coarse = isdf_data.coeff_seper;
  R_rot_coarse = isdf_data.R_rot_coarse;

  if strcmp(interp_scheme, "coarse")
    if isempty(R_rot_coarse) || size(R_rot_coarse, 1) ~= double(Nisdf)
      error('isdf:isdf_validation:BadRrot', ...
        'R_rot_coarse missing or wrong size; expected %d rows.', double(Nisdf));
    end
  end

  % Phase must use the same R_sampling_RLU as gen_tildeVq (row-aligned with coeff_seper).
  R_sampling_RLU = isdf_data.R_sampling_RLU;
  if size(R_sampling_RLU, 1) ~= double(Nisdf)
    error('isdf:isdf_validation:BadRsampling', ...
      'R_sampling_RLU rows (%d) must match coeff_seper rows (%d).', ...
      size(R_sampling_RLU, 1), double(Nisdf));
  end

  [Esum2, EsumISDF2, DiffEsum2] = isdf.validation.isdf_validate_HF(id, outDir);
  % Text summary is written by isdf.report.gen_report (called inside isdf_validate_HF) to
  %   <outDir or pwd>/isdf_validate_HF_id<id>.txt
end
