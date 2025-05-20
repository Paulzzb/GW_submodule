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
        if (options.ISDFCauchy.isISDF == true)
          Eres = gw_fullfreq_cd_res_ISDF(GWinfor, options);
          Eint = gw_fullfreq_cd_int_ISDF(GWinfor, options);
        else
          Eres = gw_fullfreq_cd_res(GWinfor, options);
          Eint = gw_fullfreq_cd_int(GWinfor, options);
        end


      else
        error('Error Not support now!');
      end
    otherwise
      error('Will never supported!!!');
  end

  % fout
  fout = options.GWCal.fileName;
  fout_flag = options.GWCal.freq_dep;
  if (fout_flag == 1); fout_flag = 0; end
  switch fout_flag 
    case 0
    %  n = options_nv - options.nv_ener+1:options_nv+options.nc_ener;
    %  ev  = GWinfor.ev * ry2ev;
    %  ev = ev(nv-nv_ener+1:nv+nc_ener);
    %  Vxc = Vxc(nv-nv_ener+1:nv+nc_ener);
      disp('We will do it later!');
    case 2
      nameConstants = fieldnames(options.Constant);
      for i = 1:numel(nameConstants)
          fieldname = nameConstants{i};
          value = options.Constant.(fieldname);    
          if ~isempty(value)
            strEval = sprintf('%s = %.16f;', fieldname, value);
            eval(strEval);
          end
      end
      % Vxc = GWinfor.Vxc * ry2ev;
      n = nv-nv_ener+1:nv+nc_ener;
      Vxc = GWinfor.Vxc;
      ev  = GWinfor.ev * ry2ev;
      ev = ev(nv-nv_ener+1:nv+nc_ener);
      Vxc = Vxc(nv-nv_ener+1:nv+nc_ener);
      Sigma = ev - Vxc + Ex + Eres + Eint;
      Eqp = ev - Vxc + Sigma;
      GWfout(2, fout, n, ev, Ex, Eres, Eint, Sigma, Vxc, Eqp);
      % save(options.GWCal.fileName, 'ev', 'Ex', 'Eres', 'Eint', 'Sigma',...
      %      'Vxc', 'Eqp');
  end


  timeGW = toc(startGW);
  fprintf('GW calculation time = %f\n', timeGW);
  end % main function

