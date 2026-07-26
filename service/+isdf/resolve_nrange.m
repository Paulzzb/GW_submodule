function [nrange1, nrange2, Nn1, Nn2] = resolve_nrange(desc, cfg_system)
%RESOLVE_NRANGE  Band ranges for vc/vn/nn from config.SYSTEM (no pool id).

  if nargin < 2 || isempty(cfg_system) || ~isstruct(cfg_system)
    error('isdf:resolve_nrange:CfgRequired', ...
      'config.SYSTEM is required (second argument).');
  end
  required_fields = {'number_bands_in_summation', ...
                     'energy_band_index_min', ...
                     'energy_band_index_max'};
  for ifd = 1:numel(required_fields)
    if ~isfield(cfg_system, required_fields{ifd})
      error('isdf:resolve_nrange:CfgField', ...
        'Missing config.SYSTEM.%s.', required_fields{ifd});
    end
  end

  desc = lower(strtrim(char(string(desc))));
  tok = regexp(desc, '^(vc|vn|nn)', 'tokens', 'once');
  if ~isempty(tok)
    desc = tok{1};
  end

  nb = double(wave_functions.get().nb);

  nocc_max = 0;
  system_data = system.get();
  nkibz = lattice.manager('k', 'get').nibz;
  for ispin = 1:system_data.nspin
    for ikibz = 1:nkibz
      idx_last = find(system_data.f(:, ikibz, ispin) > 1e-5, 1, 'last');
      if ~isempty(idx_last)
        nocc_max = max(nocc_max, idx_last);
      end
    end
  end

  nsum = min(max(1, round(double(cfg_system.number_bands_in_summation))), nb);
  nbmin = min(max(1, round(double(cfg_system.energy_band_index_min))), nb);
  nbmax = min(max(1, round(double(cfg_system.energy_band_index_max))), nb);
  if nbmin > nbmax
    error('isdf:resolve_nrange:CfgRange', ...
      'Invalid config.SYSTEM band range: min(%d) > max(%d).', ...
      int32(nbmin), int32(nbmax));
  end

  switch desc
    case 'vc'
      nrange1 = 1:nocc_max;
      nrange2 = (nocc_max + 1):nsum;
    case 'vn'
      nrange1 = 1:nocc_max;
      nrange2 = nbmin:nbmax;
    case 'nn'
      nrange1 = 1:nsum;
      nrange2 = nbmin:nbmax;
    otherwise
      error('isdf:resolve_nrange:Type', 'Unknown isdf desc ''%s''.', desc);
  end

  if isempty(nrange1) || isempty(nrange2)
    error('isdf:resolve_nrange:Range', ...
      'Empty nrange generated for desc ''%s''.', desc);
  end

  Nn1 = int32(numel(nrange1));
  Nn2 = int32(numel(nrange2));
  nrange1 = int32(nrange1);
  nrange2 = int32(nrange2);
end
