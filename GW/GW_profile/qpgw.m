function GWenergy = qpgw(GWinfo, config)

% Push functions to stack in QPlog
cleanup = QPlog_push('QP-GW');

msg = 'Main GW calculation started.';
startGW = tic;
QPlog(msg, 0);

% Initialize energy structure
GWenergy = QPenergy(GWinfo, config);
QPlog('Quasiparticle energy structure initialized.', 2);

% Step 1: Exchange part
QPlog('Calculating exchange term (Ex)...', 1);
if (config.CONTROL.enable_k_points > 1)
  GWenergy.Ex = gw_x_k(GWinfo, config);
else
  GWenergy.Ex = gw_x(GWinfo, config);
end
QPlog('Exchange term completed.', 2);

% Step 2: Frequency-dependent part
switch config.FREQUENCY.frequency_dependence 
  case -1
    QPlog('Using Hartree-Fock approximation only.', 1);
    QPlog('COHSEX calculation completed.', 2);
  case 0
    QPlog('Using COHSEX approximation (static).', 1);
    [GWenergy.Esex_x, GWenergy.Ecoh] = gw_cohsex(GWinfo, config);
    QPlog('COHSEX calculation completed.', 2);
  case 1
    QPlog('GPP method selected, but not yet implemented.', 0);
    QPerror('GPP method is under developing. Please use COHSEX or full-frequency approach.');
    % gw_gpp_kssolv(GWinfo, options);
  case 2
    msg = 'Full-frequency method selected. Generating frequency grid...';
    QPlog(msg, 1);
    config = generate_frequency(GWinfo, config);
    QPlog('Frequency grid generated.', 2);
    
    if config.FREQUENCY.frequency_dependence_method == 0 % Current Disabled.
      msg = sprintf('Real-axis approximation selected.');
      QPlog(msg, 0);
      msg = sprintf('Real-axis GW calculation is under development and currently not available.');
      QPerror(msg);
      gw_fullfreq_ra(GWinfo, config);
    elseif config.FREQUENCY.frequency_dependence_method == 2
      QPlog('Performing full-frequency GW (contour deformation method)...', 1);
      GWenergy.Eres = gw_fullfreq_cd_res(GWinfo, config);
      GWenergy.Eint = gw_fullfreq_cd_int(GWinfo, config);
    else
      msg = sprintf('config.FREQUENCY.frequency_dependence_method = %d is not supported', ...
      config.FREQUENCY.frequency_dependence_method);
      QPerror(msg);
    end
    QPlog('Full-frequency GW completed.', 2);
    % QPerror(['GW', 'fullfrequency calculation is not implemented yet.']);
  otherwise
    msg = sprintf('config.FREQUENCY.frequency_dependence %d is not supported.', ...
                   config.FREQUENCY.frequency_dependence_method);
    QPerror(msg);
end 


msg = sprintf('Total time: %.1f seconds.', toc(startGW));
QPlog(msg, 0);
msg = sprintf('GW calculation completed.');
QPlog(msg, 0);

end % EOF