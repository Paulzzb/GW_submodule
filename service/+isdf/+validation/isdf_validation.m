% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [Esum2, EsumISDF2, DiffEsum2] = isdf_validation(id)
% Run ISDF energy validation for pool id.
%
% Requires tildeVq filled (gen_tildeVq). HF text report goes to
% filename_map().isdf_report_dir via isdf.report.hf.
%
% Optional id: omit to use isdf.isdf_current().

  if nargin < 1 || isempty(id)
    id = isdf.isdf_current();
  end

  isdf_data = isdf.get(id);

  switch isdf_data.interp_scheme
    case "coarse"
      interp_scheme = "coarse";
    case "adaptive"
      interp_scheme = "adaptive";
    case "qrcp"
      interp_scheme = "qrcp";
    case "kmeans"
      interp_scheme = "kmeans";
    otherwise
      error('isdf:isdf_validation:InvalidInterpScheme', ...
        'Invalid interp_scheme: %s', isdf_data.interp_scheme);
  end

  if isempty(isdf_data.tildeVq)
    error('isdf:isdf_validation:NoTildeVq', ...
      'tildeVq is empty; run isdf.gen_tildeVq(id) first. id=%d', id);
  end

  Nisdf = isdf_data.nisdf;
  R_rot_coarse = isdf_data.R_rot_coarse;

  if strcmp(interp_scheme, "coarse")
    if isempty(R_rot_coarse) || size(R_rot_coarse, 1) ~= double(Nisdf)
      error('isdf:isdf_validation:BadRrot', ...
        'R_rot_coarse missing or wrong size; expected %d rows.', double(Nisdf));
    end
  end

  R_sampling_RLU = isdf_data.R_sampling_RLU;
  if size(R_sampling_RLU, 1) ~= double(Nisdf)
    error('isdf:isdf_validation:BadRsampling', ...
      'R_sampling_RLU rows (%d) must match coeff_seper rows (%d).', ...
      size(R_sampling_RLU, 1), double(Nisdf));
  end

  [Esum2, EsumISDF2, DiffEsum2] = isdf.validation.isdf_validate_HF(id);
end
