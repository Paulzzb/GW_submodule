function gwCalculation(GWinfor, options)
% function gwCalculation(mol, options_in)
  

  % options = GWOptions(mol, options_in);
  % GWinfor = GWinfo(mol, options.Groundstate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the jobs

  startGW = tic;
  switch options.GWCal.freq_dep
    case 0
      gw_cohsex(GWinfor, options);
    case 1
      gw_gpp_kssolv(GWinfor, options);
    case 2
      if options.GWCal.freq_dep_method == 0 % Current Disabled.
        gw_fullfreq_ra(GWinfor, options);
      elseif options.GWCal.freq_dep_method == 2
        % gw_fullfreq_cd(GWinfor, GWOptions);
        Ex = gw_x(GWinfor, options);
        Eres = gw_fullfreq_cd_res(GWinfor, options);
        Eint = gw_fullfreq_cd_int(GWinfor, options);
        Vxc = GWinfor.Vxc * options.Constant.ry2ev;
        ev  = GWinfor.ev  * options.Constant.ry2ev;
        Sigma = ev - Vxc + Ex + Eres + Eint;
        Eqp = ev - Vxc + Sigma;
        save(options.GWCal.fileName, 'ev', 'Ex', 'Eres', 'Eint', 'Sigma',...
             'Vxc', 'Eqp');
      else
        error('Error Not support now!');
      end
    otherwise
      error('Will never supported!!!');
  end
  timeGW = toc(startGW);
  fprintf('GW calculation time = %f\n', timeGW);
  end % main function

