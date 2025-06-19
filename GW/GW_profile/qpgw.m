function GWenergy = qpgw(GWinfo, config)

msg = '[QP-GW] Main GW calculation started.';
startGW = tic;
GWlog(msg, 0);

% Initialize energy structure
GWenergy = QPenergy(GWinfo, config);
GWlog('[QP-GW] Quasiparticle energy structure initialized.', 2);

% Step 1: Exchange part
GWlog('[QP-GW] Calculating exchange term (Ex)...', 1);
GWenergy.Ex = gw_x(GWinfo, config);
GWlog('[QP-GW] Exchange term completed.', 2);

% Step 2: Frequency-dependent part
switch config.FREQUENCY.frequency_dependence 
  case 0
    GWlog('[QP-GW] Using COHSEX approximation (static).', 1);
    [GWenergy.Esex_x, GWenergy.Ecoh] = gw_cohsex(GWinfo, config);
    GWlog('[QP-GW] COHSEX calculation completed.', 2);
  case 1
    GWlog('[QP-GW] GPP method selected, but not yet implemented.', 0);
    GWerror('[QP-GW] GPP method is under developing. Please use COHSEX or full-frequency approach.');
    % gw_gpp_kssolv(GWinfo, options);
  case 2
    msg = '[QP-GW] Full-frequency method selected. Generating frequency grid...';
    GWlog(msg, 1);
    config = generate_frequency(GWinfo, config);
    GWlog('[QP-GW] Frequency grid generated.', 2);
    
    if config.FREQUENCY.frequency_dependence_method == 0 % Current Disabled.
      msg = sprintf('[QP-GW] Real-axis approximation selected.');
      GWlog(msg, 0);
      msg = sprintf('[QP-GW] Real-axis GW calculation is under development and currently not available.');
      GWerror(msg);
      gw_fullfreq_ra(GWinfo, config);
    elseif config.FREQUENCY.frequency_dependence_method == 2
      GWlog('[QP-GW] Performing full-frequency GW (contour deformation method)...', 1);
      GWenergy.Eres = gw_fullfreq_cd_res(GWinfo, config);
      GWenergy.Eint = gw_fullfreq_cd_int(GWinfo, config);
    else
      msg = sprintf('[QP-GW] config.FREQUENCY.frequency_dependence_method = %d is not supported', ...
      config.FREQUENCY.frequency_dependence_method);
      GWerror(msg);
    end
    GWlog('[QP-GW] Full-frequency GW completed.', 2);
    % GWerror(['GW', 'fullfrequency calculation is not implemented yet.']);
  otherwise
    msg = sprintf('[QP-GW] config.FREQUENCY.frequency_dependence %d is not supported.', ...
                   config.FREQUENCY.frequency_dependence_method);
    GWerror(msg);
end 


msg = sprintf('[QP-GW] Total time: %.1f seconds\n', toc(startGW));
GWlog(msg, 0);
msg = sprintf('[QP-GW] GW calculation completed.\n\n');
GWlog(msg, 0)

end % EOF