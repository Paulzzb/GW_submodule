function display_input_summary(~, ~, config)

disp("========== GW Input Summary ==========")
fprintf(" Groundstate Type     : %s\n", config.CONTROL.groundstate_type);
fprintf(" Groundstate Dir      : %s\n", config.CONTROL.groundstate_dir);
fprintf(" Storage Dir          : %s\n", config.CONTROL.storage_dir);
fprintf(" Output Dir           : %s\n", config.CONTROL.output_dir);
fprintf(" Output File          : %s\n", config.CONTROL.outfile);
fprintf(" Prefix               : %s\n", config.CONTROL.prefix);

disp("----------- ISDF Settings -----------")
fprintf(" ISDF Enabled         : %d\n", config.ISDF.isisdf);
fprintf(" ISDF Ratio           : %.2f\n", config.ISDF.isdf_ratio);
fprintf(" ISDF Type1 Ratio     : %.2f\n", config.ISDF.isdf_ratio_type1);
fprintf(" ISDF Type2 Ratio     : %.2f\n", config.ISDF.isdf_ratio_type2);
fprintf(" ISDF Type3 Ratio     : %.2f\n", config.ISDF.isdf_ratio_type3);
fprintf(" Adaptive Thres T1/T2/T3: %.2e / %.2e / %.2e\n", ...
  config.ISDF.adaptive_threshold_type1, config.ISDF.adaptive_threshold_type2, ...
  config.ISDF.adaptive_threshold_type3);
fprintf(" Adaptive NumAdd T1/T2/T3: %d / %d / %d\n", ...
  int32(config.ISDF.adaptive_num_add_type1), int32(config.ISDF.adaptive_num_add_type2), ...
  int32(config.ISDF.adaptive_num_add_type3));
fprintf(" Adaptive CandR T1/T2/T3: %.2f / %.2f / %.2f\n", ...
  config.ISDF.adaptive_candidate_ratio_type1, config.ISDF.adaptive_candidate_ratio_type2, ...
  config.ISDF.adaptive_candidate_ratio_type3);
fprintf(" Adaptive MaxFrac T1/T2/T3: %.2f / %.2f / %.2f\n", ...
  config.ISDF.adaptive_max_add_frac_type1, config.ISDF.adaptive_max_add_frac_type2, ...
  config.ISDF.adaptive_max_add_frac_type3);
fprintf(" Adaptive MaxCond T1/T2/T3: %.2e / %.2e / %.2e\n", ...
  config.ISDF.adaptive_max_cond_type1, config.ISDF.adaptive_max_cond_type2, ...
  config.ISDF.adaptive_max_cond_type3);
if isfield(config.ISDF, 'adaptive_use_cond_guard')
  fprintf(" Adaptive CondGuard    : %d\n", logical(config.ISDF.adaptive_use_cond_guard));
end
if isfield(config.ISDF, 'validate_hf')
  fprintf(" ISDF validate_hf      : %d\n", logical(config.ISDF.validate_hf));
end
if isfield(config.ISDF, 'inv_strategy')
  fprintf(" ISDF inv_strategy     : %s\n", char(string(config.ISDF.inv_strategy)));
end
if isfield(config.ISDF, 'inv_param')
  fprintf(" ISDF inv_param        : %.2e\n", double(config.ISDF.inv_param));
end
if isfield(config.ISDF, 'inv_ratio')
  fprintf(" ISDF inv_ratio        : %.4f\n", double(config.ISDF.inv_ratio));
end
if isfield(config.ISDF, 'auto_inv_param')
  fprintf(" ISDF auto_inv_param   : %d\n", logical(config.ISDF.auto_inv_param));
end
if isfield(config.ISDF, 'order')
  fprintf(" ISDF order (tildeVq)  : %d\n", int32(config.ISDF.order));
end
if isfield(config.ISDF, 'chol_maxit')
  fprintf(" ISDF chol_maxit        : %d\n", int32(config.ISDF.chol_maxit));
end

disp("----------- CUTOFFS Settings -----------")
fprintf(" COULOMB TRUNCATION   : %3d\n", config.CUTOFFS.coulomb_truncation_method);
fprintf(" TRUNCATION RADIUS    : %.2f (Ry)\n", config.CUTOFFS.coulomb_truncation_parameter);
fprintf(" COULOMB CUTOFF RADIUS: %.2f (Ry)\n", config.CUTOFFS.coulomb_cutoff);
if config.FREQUENCY.frequency_dependence == 1
  fprintf(" DENSITY CUTOFF RADIUS: %.2f (Ry)\n", config.CUTOFFS.density_cutoff);
end

disp("-------- Frequency Settings ---------")
fprintf(" Frequency Dependence : %d\n", config.FREQUENCY.frequency_dependence);
if config.FREQUENCY.frequency_dependence == 2
  fprintf(" Frequency Dependence : %d\n", config.FREQUENCY.frequency_dependence);
  % fprintf(" Delta Frequency      : %.2f\n", config.FREQUENCY.delta_frequency);
  % fprintf(" Eta                  : %.1e\n", config.FREQUENCY.eta);
  fprintf("MORE INFORMATION NEEDED TO DISPLAY")
end

disp("----------- COHSEX Settings -----------")
if isfield(config.COHSEX, 'exact_ch')
  fprintf(" COHSEX exact_ch      : %d\n", logical(config.COHSEX.exact_ch));
end
if isfield(config.COHSEX, 'ex_use_which_isdf')
  fprintf(" COHSEX Ex ISDF slot  : %s\n", char(string(config.COHSEX.ex_use_which_isdf)));
end

if isfield(config, 'SUPERCELL')
  disp("----------- SUPERCELL Settings -----------")
  if isfield(config.SUPERCELL, 'is_supercell')
    fprintf(" SUPERCELL enabled    : %d\n", logical(config.SUPERCELL.is_supercell));
  end
  if isfield(config.SUPERCELL, 'use_sc_isdf')
    fprintf(" SUPERCELL use SC_ISDF: %d\n", logical(config.SUPERCELL.use_sc_isdf));
  end
  if isfield(config.SUPERCELL, 'sc_adaptive')
    fprintf(" SUPERCELL sc_adaptive : %d\n", logical(config.SUPERCELL.sc_adaptive));
  end
  if isfield(config.SUPERCELL, 'k1')
    fprintf(" SUPERCELL k          : [%d %d %d]\n", ...
      int32(config.SUPERCELL.k1), int32(config.SUPERCELL.k2), int32(config.SUPERCELL.k3));
  end
  if isfield(config.SUPERCELL, 'isdf_source_dir') && ~isempty(config.SUPERCELL.isdf_source_dir)
    fprintf(" SUPERCELL ISDF source: %s\n", char(string(config.SUPERCELL.isdf_source_dir)));
  end
end

if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'groundstate_type') ...
    && strcmpi(char(string(config.CONTROL.groundstate_type)), 'formal')
  disp("----------- FORMAL benchmark Settings -----------")
  fprintf(" groundstate_type   : formal (synthetic wf/Eo, nsym=1)\n");
  if isfield(config, 'FORMAL')
    if isfield(config.FORMAL, 'wf_seed')
      fprintf(" FORMAL wf_seed      : %d\n", round(config.FORMAL.wf_seed));
    end
    if isfield(config.FORMAL, 'eo_e0')
      fprintf(" FORMAL eo_e0        : %.6f Ry\n", config.FORMAL.eo_e0);
    end
    if isfield(config.FORMAL, 'eo_delta')
      fprintf(" FORMAL eo_delta     : %.6f Ry\n", config.FORMAL.eo_delta);
    end
  end
end

disp("====================================")
end % EOF
