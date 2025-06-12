function qp_driver(input_dir)
% qp_driver -> driver to perform quasi-particle calculation

dir = config.CONTROL.storage_dir;
def = filename_map();
fName = fullfile(dir, def.GWinput);
TEMP = load(fName);


GWoptions = construct_GWOptions(data, config);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the jobs

  startGW = tic;
  GWenergy = QPenergy(GWinfor, options);
  
  % Main calculation of GW
  GWenergy.Ex = gw_x(GWinfor, config, options.ISDFCauchy);
  switch options.GWCal.freq_dep
    case 0
      [GWenergy.Esex_x, GWenergy.Ecoh] = gw_cohsex(GWinfor, options);
    case 1
      gw_gpp_kssolv(GWinfor, options);
    case 2
      if options.GWCal.freq_dep_method == 0 % Current Disabled.
        gw_fullfreq_ra(GWinfor, options);
      elseif options.GWCal.freq_dep_method == 2
        % gw_fullfreq_cd(GWinfor, GWOptions);
        % if (options.ISDFCauchy.isISDF == true)
        %   GWenergy.Eres = gw_fullfreq_cd_res_ISDF(GWinfor, options);
        %   GWenergy.Eint = gw_fullfreq_cd_int_ISDF(GWinfor, options);
        % else
        %   GWenergy.Eres = gw_fullfreq_cd_res(GWinfor, options);
        %   GWenergy.Eint = gw_fullfreq_cd_int(GWinfor, options);
        % end
        GWenergy.Eres = gw_fullfreq_cd_res2(GWinfor, options);
        GWenergy.Eint = gw_fullfreq_cd_int2(GWinfor, options);
      else
        error('Error Not support now!');
      end
    otherwise
      error('Will never supported!!!');
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

