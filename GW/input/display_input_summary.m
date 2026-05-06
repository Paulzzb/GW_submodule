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
fprintf(" Adaptive Thres T1/T2 : %.2e / %.2e\n", ...
  config.ISDF.adaptive_threshold_type1, config.ISDF.adaptive_threshold_type2);
fprintf(" Adaptive NumAdd T1/T2: %d / %d\n", ...
  int32(config.ISDF.adaptive_num_add_type1), int32(config.ISDF.adaptive_num_add_type2));
fprintf(" Adaptive CandR T1/T2 : %.2f / %.2f\n", ...
  config.ISDF.adaptive_candidate_ratio_type1, config.ISDF.adaptive_candidate_ratio_type2);
fprintf(" Adaptive MaxFrac T1/T2: %.2f / %.2f\n", ...
  config.ISDF.adaptive_max_add_frac_type1, config.ISDF.adaptive_max_add_frac_type2);

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

disp("====================================")
end % EOF