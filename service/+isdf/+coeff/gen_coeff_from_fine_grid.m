% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11

function gen_coeff_from_fine_grid(id)
%GEN_COEFF_FROM_FINE_GRID  Fill coeff_seper from existing fine-grid sampling indices.
%
% Expects isdf_data.bundle_struct.fine_grid_lin (or R_sampling_RLU rows) after
% SC_ISDF replication. Does not change the sampling grid itself.

  isdf_data = isdf.get(id);
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');

  if isfield(isdf_data, 'bundle_struct') && isfield(isdf_data.bundle_struct, 'fine_grid_lin') ...
      && ~isempty(isdf_data.bundle_struct.fine_grid_lin)
    lin = int32(isdf_data.bundle_struct.fine_grid_lin(:));
  else
    lin = isdf.SC_ISDF_r_sampling_to_lin(isdf_data);
  end

  Nmu = numel(lin);
  nb = wf_data.nb;
  nbz = k_data.nbz;
  nspin = wf_data.nspin;
  coeff = zeros(double(Nmu), nb, nbz, nspin);

  ispin = 1;
  for ikbz = 1:nbz
    ikibz = k_data.bz2ibz(ikbz, 1);
    ikrot = k_data.bz2rot(ikbz, 1);
    for ib = 1:nb
      isc = [ib, ikibz, ikrot, ispin];
      wf = wave_functions.WF_apply_symm(isc);
      coeff(:, ib, ikbz, ispin) = wf(lin);
    end
  end

  isdf_data.coeff_seper = coeff;
  isdf_data.nisdf = int32(Nmu);
  if isfield(isdf_data, 'bundle_struct') && isstruct(isdf_data.bundle_struct)
    isdf_data.bundle_struct.N_sampling = int32(Nmu);
    isdf_data.bundle_struct.N_bundle = int32(Nmu);
    isdf_data.bundle_struct.sampling2bundle = int32((1:Nmu).');
    isdf_data.bundle_struct.R_grid_bundle = isdf_data.R_sampling_RLU(1:Nmu, :);
    isdf_data.bundle_struct.WF_bundle = coeff;
    if ~isfield(isdf_data.bundle_struct, 'N_coarse') || isempty(isdf_data.bundle_struct.N_coarse)
      isdf_data.bundle_struct.N_coarse = isdf_data.N_coarse;
    end
  end
  isdf.save2mod(isdf_data, id);
end
