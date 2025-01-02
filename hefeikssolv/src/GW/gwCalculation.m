function gwCalculation(mol, options_in)
  

  options = GWOptions(mol, options_in);
  ksinfor = ksinfo(mol, options.Groundstate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the jobs

  startGW = tic;
  switch options.GWCal.freq_dep
    case 0
      gw_cohsex(ksinfor, options);
    case 1
      gw_gpp_kssolv(ksinfor, options);
    case 2
      if options.frequency_dependence_method == 0 
        gw_fullfreq_ra(ksinfor, options);
      elseif options.frequency_dependence_method == 2
        % gw_fullfreq_cd(ksinfor, GWOptions);
        Ex = gw_x(ksinfor, GWOptions);
        Eres = gw_fullfreq_cd_res(ksinfor, GWOptions);
        Eint = gw_fullfreq_cd_int(ksinfor, GWOptions);
        Vxc = ksinfo.Vxc * options.Constant.ry2ev;
        ev  = ksinfo.ev  * options.Constant.ry2ev;
        Sigma = ev - Vxc + Ex + Eres + Eint
      else
        error('Error Not support now!');
      end
    otherwise
      error('Will never supported!!!');
  end
  timeGW = toc(startGW);
  fprintf('GW calculation time = %f\n', timeGW);
  end % main function

