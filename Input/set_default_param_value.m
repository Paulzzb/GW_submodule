function config = set_default_param_value(config, data)
% This function fills in default values for missing config
% parameters with the help of system information
% Call after read_input_param

defaults = default_param_values();

blocks = fieldnames(defaults);
for i = 1:numel(blocks)
    blk = blocks{i};
    if ~isfield(config, blk)
        config.(blk) = struct();
    end
    keys = fieldnames(defaults.(blk));
    for j = 1:numel(keys)
        key = keys{j};
        if ~isfield(config.(blk), key)
            config.(blk).(key) = defaults.(blk).(key);
        end
    end
end

% Fill in default values for missing parameters depending on data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters
if (config.SYSTEM.number_bands_in_summation < 0)
  config.SYSTEM.number_bands_in_summation = length(data.ev);
end

if (config.SYSTEM.energy_band_index_min < 0)
  config.SYSTEM.energy_band_index_min= 1;
end

if (config.SYSTEM.energy_band_index_max < 0)
  config.SYSTEM.energy_band_index_max = length(data.ev);
end

if (config.CUTOFFS.coulomb_cutoff < 0)
  config.CUTOFFS.coulomb_cutoff = data.reciprocal_grid_info.wfncut;
end

if (config.CUTOFFS.coulomb_cutoff > data.reciprocal_grid_info.wfncut)
  error("&SYSTEM->coulomb_cutoff = %8.2f Ry should be no bigger than wavefuction cutoff = %8.2f Ry.", ...
        config.CUTOFFS.coulomb_cutoff, data.reciprocal_grid_info.wfncut);
end

if (config.CUTOFFS.density_cutoff < 0)
  config.CUTOFFS.density_cutoff = config.CUTOFFS.density_cutoff * 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters with respect to the GW approximating strategy
if (config.FREQUENCY.frequency_dependence == 0)
  % COHSEX
  ;
elseif (config.FREQUENCY.frequency_dependence == 1)
  % GPP, complete it later
  ;
elseif (config.FREQUENCY.frequency_dependence == 2)
  % FULL FREQUENCY
  if (config.FREQUENCY.frequency_low_cutoff < 0)
    % Now only insulaters are supported, and consider spin degeneracy only.
    nv = data.sys.ne / 2;
    ev = data.ev;
    nbmax = max(config.SYSTEM.energy_band_index_max, ...
                config.SYSTEM.number_bands_in_summation);
    tmp1 = ev(nv) - ev(1); tmp2 = ev(nbmax) - ev(nv+1);

    config.FREQUENCY.frequency_low_cutoff = max(tmp1, tmp2); 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters with respect to the ISDF 
if (config.ISDF.isisdf > 0)
  if (config.ISDF.isdf_ratio < 0)
    error('&ISDF->isdf_ratio should be set to a positive value!');
  end
  if (config.ISDF.isdf_ratio_type1 < 0)
    config.ISDF.isdf_ratio_type1 = config.ISDF.isdf_ratio;
  end
  if (config.ISDF.isdf_ratio_type2 < 0)
    config.ISDF.isdf_ratio_type2 = config.ISDF.isdf_ratio;
  end
  if (config.ISDF.isdf_ratio_type3 < 0)
    config.ISDF.isdf_ratio_type3 = config.ISDF.isdf_ratio;
  end
  if isempty(config.ISDF.sys)
    config.ISDF.sys = data.sys;
  end
end


end