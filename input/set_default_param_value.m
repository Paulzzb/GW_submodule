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

% if (config.CUTOFFS.coulomb_cutoff > data.reciprocal_grid_info.wfncut)
%   error("&SYSTEM->coulomb_cutoff = %8.2f Ry should be no bigger than wavefuction cutoff = %8.2f Ry.", ...
%         config.CUTOFFS.coulomb_cutoff, data.reciprocal_grid_info.wfncut);
% end

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
    ha2ry = 2.0;
    ev = ha2ry * data.ev; % Remember input ev is Ha, while we hope all parameters are in Ry.
    nbmax = max(config.SYSTEM.energy_band_index_max, ...
                config.SYSTEM.number_bands_in_summation);
    tmp1 = ev(nv) - ev(1); tmp2 = ev(nbmax) - ev(nv+1);

    config.FREQUENCY.frequency_low_cutoff = max(tmp1, tmp2) * 1.25; 
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
  if (config.ISDF.adaptive_threshold_type1 < 0)
    config.ISDF.adaptive_threshold_type1 = 2e-4;
  end
  if (config.ISDF.adaptive_threshold_type2 < 0)
    config.ISDF.adaptive_threshold_type2 = 2e-4;
  end
  if (config.ISDF.adaptive_threshold_type3 < 0)
    config.ISDF.adaptive_threshold_type3 = 2e-4;
  end
  if (config.ISDF.adaptive_num_add_type1 < 0)
    config.ISDF.adaptive_num_add_type1 = 16;
  end
  if (config.ISDF.adaptive_num_add_type2 < 0)
    config.ISDF.adaptive_num_add_type2 = 16;
  end
  if (config.ISDF.adaptive_num_add_type3 < 0)
    config.ISDF.adaptive_num_add_type3 = 16;
  end
  if (config.ISDF.adaptive_candidate_ratio_type1 < 0)
    config.ISDF.adaptive_candidate_ratio_type1 = 2.0;
  end
  if (config.ISDF.adaptive_candidate_ratio_type2 < 0)
    config.ISDF.adaptive_candidate_ratio_type2 = 2.0;
  end
  if (config.ISDF.adaptive_candidate_ratio_type3 < 0)
    config.ISDF.adaptive_candidate_ratio_type3 = 2.0;
  end
  if (config.ISDF.adaptive_max_add_frac_type1 < 0)
    config.ISDF.adaptive_max_add_frac_type1 = 1.0;
  end
  if (config.ISDF.adaptive_max_add_frac_type2 < 0)
    config.ISDF.adaptive_max_add_frac_type2 = 1.0;
  end
  if (config.ISDF.adaptive_max_add_frac_type3 < 0)
    config.ISDF.adaptive_max_add_frac_type3 = 1.0;
  end
  if (config.ISDF.adaptive_max_cond_type1 < 0)
    config.ISDF.adaptive_max_cond_type1 = 1e6;
  end
  if (config.ISDF.adaptive_max_cond_type2 < 0)
    config.ISDF.adaptive_max_cond_type2 = 1e6;
  end
  if (config.ISDF.adaptive_max_cond_type3 < 0)
    config.ISDF.adaptive_max_cond_type3 = 1e6;
  end
  if ~(isfinite(config.ISDF.inv_ratio) && config.ISDF.inv_ratio >= 0 && config.ISDF.inv_ratio <= 1)
    config.ISDF.inv_ratio = 0.5;
  end
  if isempty(config.ISDF.sys)
    config.ISDF.sys = data.sys;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supercell / SC_ISDF
if ~isfield(config, 'SUPERCELL')
  config.SUPERCELL = default_param_values().SUPERCELL;
end
if config.SUPERCELL.k1 < 1
  config.SUPERCELL.k1 = 1;
end
if config.SUPERCELL.k2 < 1
  config.SUPERCELL.k2 = 1;
end
if config.SUPERCELL.k3 < 1
  config.SUPERCELL.k3 = 1;
end
if config.SUPERCELL.use_sc_isdf && ~config.SUPERCELL.is_supercell
  config.SUPERCELL.is_supercell = true;
end
if config.SUPERCELL.use_sc_isdf
  src = strtrim(char(string(config.SUPERCELL.isdf_source_dir)));
  if isempty(src)
    error('SUPERCELL.use_sc_isdf=true requires SUPERCELL.isdf_source_dir.');
  end
  if ~isfolder(src)
    error('SUPERCELL.isdf_source_dir does not exist: %s', src);
  end
end
if config.SUPERCELL.sc_adaptive && ~config.SUPERCELL.use_sc_isdf
  error('SUPERCELL.sc_adaptive=true requires SUPERCELL.use_sc_isdf=true.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formal benchmark (synthetic groundstate + Gamma COHSEX path)
if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'groundstate_type') ...
    && strcmpi(char(string(config.CONTROL.groundstate_type)), 'formal')
  if config.FREQUENCY.frequency_dependence ~= -2
    error('groundstate_type=''formal'' requires FREQUENCY.frequency_dependence = -2.');
  end
  if ~config.SUPERCELL.use_sc_isdf
    error('groundstate_type=''formal'' requires SUPERCELL.use_sc_isdf = true.');
  end
  config.SUPERCELL.is_supercell = true;
  config.CONTROL.enable_k_points = false;
  config.ISDF.validate_hf = false;
  if ~isfield(config, 'COHSEX') || ~isstruct(config.COHSEX)
    config.COHSEX = default_param_values().COHSEX;
  end
  config.COHSEX.exact_ch = true;
  config.ISDF.compute_vc = true;
  config.ISDF.compute_vn = true;
end

end
