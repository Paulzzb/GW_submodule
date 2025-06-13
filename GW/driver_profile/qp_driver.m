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
      GWerror(['GW', 'fullfrequency calculation is not implemented yet.']);
      % if options.GWCal.freq_dep_method == 0 % Current Disabled.
      %   gw_fullfreq_ra(GWinfo, options);
      % elseif options.GWCal.freq_dep_method == 2
      %   % gw_fullfreq_cd(GWinfo, GWOptions);
      %   % if (options.ISDFCauchy.isISDF == true)
      %   %   GWenergy.Eres = gw_fullfreq_cd_res_ISDF(GWinfo, options);
      %   %   GWenergy.Eint = gw_fullfreq_cd_int_ISDF(GWinfo, options);
      %   % else
      %   %   GWenergy.Eres = gw_fullfreq_cd_res(GWinfo, options);
      %   %   GWenergy.Eint = gw_fullfreq_cd_int(GWinfo, options);
      %   % end
      %   GWenergy.Eres = gw_fullfreq_cd_res2(GWinfo, options);
      %   GWenergy.Eint = gw_fullfreq_cd_int2(GWinfo, options);
      % else
      %   error('Error Not support now!');
      % end
    otherwise
      GWerror('config.FREQUENCY.frequency_dependence is not supported');
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

