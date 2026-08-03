% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/13 ZZ

function set_nrange(id, cfg_system)
  % Get nrange for "vc"/"vn"/"nn" with cache in isdf_data.
  if nargin < 2 || isempty(cfg_system) || ~isstruct(cfg_system)
    error('isdf_set_nrange:CfgRequired', ...
      'config.SYSTEM is required (second argument).');
  end
  required_fields = {'number_bands_in_summation', ...
                     'energy_band_index_min', ...
                     'energy_band_index_max'};
  for ifd = 1:numel(required_fields)
    if ~isfield(cfg_system, required_fields{ifd})
      error('set_nrange:CfgField', ...
        'Missing config.SYSTEM.%s.', required_fields{ifd});
    end
  end

  isdf_data = isdf.get(id);
  if ~isempty(isdf_data.nrange1) && ~isempty(isdf_data.nrange2)
    return;
  end

  desc = char(string(isdf_data.desc));
  [nrange1, nrange2, Nn1, Nn2] = isdf.resolve_nrange(desc, cfg_system);

  isdf_data.Nnrange1 = Nn1;
  isdf_data.Nnrange2 = Nn2;
  isdf_data.nrange1 = nrange1;
  isdf_data.nrange2 = nrange2;
  isdf.save2mod(isdf_data, id);
end
