function qp_driver(input_dir)
% qp_driver -> driver to perform quasi-particle calculation

def = filename_map();
fName = fullfile(input_dir, def.GWinput);
TEMP = load(fName);

GWinfo = TEMP.GWgroundstate;
config = TEMP.config;
GWoptions = TEMP.GWoptions;

init_log(config.CONTROL.log_level)


% Prepare GWinfo.psir
GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the jobs

  startGW = tic;
  GWenergy = QPenergy(GWinfo, config);
  
  % Main calculation of GW
  GWenergy.Ex = gw_x(GWinfo, config);
  switch config.FREQUENCY.frequency_dependence 
    case 0
      [GWenergy.Esex_x, GWenergy.Ecoh] = gw_cohsex(GWinfo, config);
    case 1
      GWerror(['GW', 'GPP calculation is not implemented yet.']);
      % gw_gpp_kssolv(GWinfo, options);
    case 2
      config = generate_frequency(GWinfo, config);
      
      if config.FREQUENCY.frequency_dependence_method == 0 % Current Disabled.
        msg = sprintf('Implementing real axis approximation');
        GWlog(msg, 0);
        msg = sprintf('GW calculation under real axis approximation is developing');
        GWerror(msg);
        gw_fullfreq_ra(GWinfo, config);
      elseif config.FREQUENCY.frequency_dependence_method == 2
        GWenergy.Eres = gw_fullfreq_cd_res(GWinfo, config);
        GWenergy.Eint = gw_fullfreq_cd_int(GWinfo, config);
      else
        msg = sprintf('config.FREQUENCY.frequency_dependence_method = %d is not supported', ...
        config.FREQUENCY.frequency_dependence_method);
        GWerror(msg);
      end
      % GWerror(['GW', 'fullfrequency calculation is not implemented yet.']);
    otherwise
      GWerror('config.FREQUENCY.frequency_dependence = %d is not supported', ...
               config.FREQUENCY.frequency_dependence);
  end

  % Degeneracy
  GWenergy = shiftenergy(GWenergy);
  
  % Calculate based on energies info only
  % % Dyn with respect to frequency !!
  GWenergy = getEqp(GWenergy);
  % Output
  GWfout(GWenergy);

  % End of the function
  timeGW = toc(startGW);
  fprintf('GW calculation time = %f\n', timeGW);
  end % main function

